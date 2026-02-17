#from __init__ import *
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
import mylib.utils.galfit as galfit
import time

from scipy.ndimage.measurements import center_of_mass
import emcee
import pickle
import corner

#dir_root = os.path.abspath(os.getcwd()+'/../')
#dir_data = dir_root + '/data'
#dir_code = dir_root + '/code'
#dir_exe = dir_root + '/code'

glafic_exe =  './glafic_mac'

def PSF_command(paras_psf1, paras_psf2, psf_ratio):
    FWHM1, e1, theta_e1, beta1 = paras_psf1
    FWHM2, e2, theta_e2, beta2 = paras_psf2

    command = 'psf %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n'\
        %(FWHM1, e1, theta_e1, beta1, FWHM2, e2, theta_e2, beta2, psf_ratio)
    return command


class OneObject(object):
    def __init__(self, redshift, mass_component=[], light_component=[]):
        self.Masslist = mass_component.copy()
        self.Lightlist = light_component.copy()
        self.redshift = redshift
        self.get_component_numbers()

    def get_component_numbers(self):
        self.n_mass = len(self.Masslist)
        self.n_light = len(self.Lightlist)

    def add_mass_component(self, mass_component):
        self.Masslist.append(mass_component)
        self.get_component_numbers()

    def reset_mass_component(self):
        self.Masslist = []
        self.get_component_numbers()

    def add_light_component(self, light_component):
        self.Lightlist.append(light_component)
        self.get_component_numbers()

    def reset_light_component(self):
        self.Lightlist = []
        self.get_component_numbers()


class MassObject(object):
    def __init__(self):
        do_something = 1

class LightObject(object):
    def __init__(self):
        do_something = 1

class LensSystem(object):

    def __init__(self, lensobj, sourceobj,\
                 kwargs_cosmo={'omega':0.3, 'lambda':0.7, 'hubble':0.7}):
        self.lens = lensobj
        self.source = sourceobj
        self.kwargs_cosmo = kwargs_cosmo
        self.image_dat = self.findimage()

    def _get_primary_command(self, **kwargs):

        # Primary Parameters #
        omega_m = self.kwargs_cosmo['omega']
        omega_l = self.kwargs_cosmo['lambda']
        hubble = self.kwargs_cosmo['hubble']

        zl = self.lens.redshift

        xmin, xmax = kwargs['x_range']
        ymin, ymax = kwargs['y_range']
        prefix = kwargs['prefix']
        pix_ext = kwargs['pix_ext']
        pix_poi = kwargs['pix_poi']

        primary_param_str =\
'''
# primary parameters
omega    %.3f
lambda   %.3f
weos     -1.000
hubble   %.3f
zl       %.3f
prefix   %s
xmin     %.3f
ymin     %.3f
xmax     %.3f
ymax     %.3f
pix_ext  %.3f
pix_poi  %.3f
maxlev   6
'''%(omega_m, omega_l, hubble, zl, prefix, xmin, ymin, xmax, ymax,\
             pix_ext, pix_poi)

        return primary_param_str

    def _get_secondary_command(self, kwargs_secondary):
        command = ''
        for item in kwargs_secondary.items():
            key, value = item
            command += '%s  %s\n'%(key, value)

        return command

    def _get_lens_command(self):

        lensstr = ''

        # lens model command #
        for idx_lens in range(self.lens.n_mass):
            lens_mass = self.lens.Masslist[idx_lens]
            lens_x = lens_mass['x_center']
            lens_y = lens_mass['y_center']
            lens_ellipticity = lens_mass['ellipticity']
            lens_pa = lens_mass['pa']
            lens_type = lens_mass['Type']

            if lens_type == 'nfw':
                key_mass = lens_mass['mass']
                key_6 = lens_mass['c']
                key_7 = 0

            elif lens_type == 'sers':
                key_mass = lens_mass['mass']
                key_6 = lens_mass['re']
                key_7 = lens_mass['n']

            elif lens_type == 'pow':
                key_mass = lens_mass['mass']
                key_6 = lens_mass['r_ein']
                key_7 = lens_mass['gamma']

            elif lens_type == 'sie':
                key_mass = lens_mass['sigma']
                key_6 = lens_mass['r_core']
                key_7 = 0

            else:
                raise UnkownMassProfileError(\
                                'Unkown Lens Mass Profile %s'%lens_type)

            lensstr += 'lens %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n'\
                    %(lens_type, key_mass, lens_x, lens_y,\
                      lens_ellipticity, lens_pa, key_6, key_7)

        return lensstr

    def _get_source_command(self, islens=False):
        # source model command #
        sourcestr = ''
        if islens:
            source_to_model = self.lens
        else:
            source_to_model = self.source

        for idx_source in range(source_to_model.n_light):
            source = source_to_model.Lightlist[idx_source]
            source_z = source_to_model.redshift
            if islens:
                source_z +=0.1

            source_flux = source['flux']
            source_x = source['x_center']
            source_y = source['y_center']
            source_type = source['Type']

            if source_type=='point':

                sourcestr += 'extend %s %.5f %.5f %.5f %.5f\
                        %.5f %.5f %.5f %.5f\n'\
                    %(source_type, source_z, source_flux, source_x, source_y,\
                     0, 0, 0, 0)

            elif source_type=='sersic':
                source_ellipticity = source['ellipticity']
                source_pa = source['pa']
                source_re = source['re']
                source_n = source['n']

                sourcestr +=\
                    'extend %s %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n'\
                    %(source_type, source_z, source_flux, source_x, source_y,\
                    source_ellipticity, source_pa, source_re, source_n)

            else:
                raise UnkownLightProfileError(\
                                'Unkown Source Light Profile %s'%lens_type)
        return sourcestr

    def _get_pointsource_command(self):
        pointsourcestr = ''

        source_to_model = self.source
        for idx_source in range(source_to_model.n_light):
            source = source_to_model.Lightlist[idx_source]
            source_z = source_to_model.redshift
            source_x = source['x_center']
            source_y = source['y_center']

            pointsourcestr += 'point %.5f %.5f %.5f\n'\
                    %(source_z, source_x, source_y)

        return pointsourcestr

    def _get_object_command(self):
        num_lens = self.lens.n_mass
        num_ps = self.source.n_light
        num_ext = 0

        object_str = 'startup %d %d %d\n'%(num_lens, num_ext, num_ps)

        lensstr = self._get_lens_command()
        sourcestr = self._get_pointsource_command()
        psfstr = ''

        object_str = object_str + lensstr + sourcestr + psfstr\
                + 'end_startup\n'

        return object_str

    def _get_optimization_command(self, command=''):
        if len(command)==0:
            command = ''
        else:
            command = 'start_setopt\n'+command+'\nend_setopt\n'
        return command

    def _get_operation_command(self, operation_list):
        operation_command = 'start_command\n'
        for operation in operation_list:
            operation_command += '%s\n'%(operation)

        operation_command += 'quit\n'
        return operation_command

    def glafic_command(self, prefix, x_range, y_range, pix_ext, pix_poi,\
                       operation_list, kwargs_secondary):

        primary_command = self._get_primary_command(\
                    pix_ext=pix_ext, prefix=prefix, pix_poi=pix_poi,\
                    x_range=x_range, y_range=y_range)
        secondary_command = self._get_secondary_command(kwargs_secondary)
        object_command = self._get_object_command()
        optimization_command = self._get_optimization_command()
        operation_command = self._get_operation_command(operation_list)

        allcommand = primary_command + secondary_command\
                + object_command +optimization_command + operation_command

        return allcommand

    def run_command(self, prefix, x_range, y_range, pix_ext, pix_poi,\
                 command_list, **kwargs):

        tmp_input_file = 'tmp%s%d.input'%(prefix, np.random.randint(100))
        allcommand = self.glafic_command(prefix, x_range, y_range,\
                                pix_ext, pix_poi, command_list, kwargs)
        f = open(dir_exe + '/' + tmp_input_file, 'w+')
        f.write(allcommand)
        f.close()

        os.system(glafic_exe + ' %s/%s > /dev/null'%(dir_exe, tmp_input_file))
#        os.system('rm %s/%s'%(dir_exe, tmp_input_file))

    def plot_lensconfig(self, output_image, x_range, y_range):
        image_dat = self.image_dat
        crit_dat = self.findcrit()

        image_x, image_y, image_mu, image_dt = image_dat
        source_x = self.source.Lightlist[0]['x_center']
        source_y = self.source.Lightlist[0]['y_center']

        xi1, yi1, xs1, ys1, xi2, yi2, xs2, ys2 = crit_dat

        fig, ax = plt.subplots(figsize=[5, 5])
        ax.plot(image_x, image_y, 'c+')
        ax.plot(source_x, source_y, 'co')
        for index in range(len(xi1)):
            ax.plot([xi1[index], xi2[index]], [yi1[index], yi2[index]],\
                    'r', lw=0.5)
            ax.plot([xs1[index], xs2[index]], [ys1[index], ys2[index]],\
                    'b', lw=0.5)

        ax.set_xlabel('x (arcsec)')
        ax.set_xlim(x_range)
        ax.set_ylabel('y (arcsec)')
        ax.set_ylim(y_range)
        ax.set_aspect('equal')
        plt.tight_layout()
        fig.savefig(output_image)

    def findimage(self, x_range=[-10, 10], y_range=[-10, 10],\
                  pix_ext=0.05, pix_poi=0.05, prefix='out', **kwargs):
        self.run_command(prefix, x_range, y_range, pix_ext, pix_poi,\
                      command_list=['findimg'], **kwargs)
        data = np.loadtxt(prefix + '_point.dat')
        return data[1:].T

    def findcrit(self, x_range=[-10, 10], y_range=[-10, 10],\
                  pix_ext=0.01, pix_poi=0.01, prefix='out', **kwargs):
        prefix = 'out'
        self.run_command(prefix, x_range, y_range,  pix_ext, pix_poi,\
                      command_list=['writecrit'], **kwargs)
        data = np.loadtxt(prefix + '_crit.dat')
        return data.T

    def update_sourcemag(self, mag, index=0):
        self.source.Lightlist[index]['Mag'] = mag

    def update_lensmag(self, mag, index=0):
        self.lens.Lightlist[index]['Mag'] = mag

    def lens_info(self, img_size, pix_scale):
        lens_x_pix_list = []
        lens_y_pix_list = []
        mag_list = []
        re_list = []
        n_list = []
        pa_list = []
        ar_list = []

        for index in range(self.lens.n_light):
            lens_x = self.lens.Lightlist[index]['x_center']
            lens_y = self.lens.Lightlist[index]['y_center']

            lens_x_pix = (lens_x + img_size/2)/pix_scale + 1
            lens_y_pix = (lens_y + img_size/2)/pix_scale + 1

            mag = self.lens.Lightlist[index]['Mag']
            re = self.lens.Lightlist[index]['re']/pix_scale
            n = self.lens.Lightlist[index]['n']
            pa = self.lens.Lightlist[index]['pa']
            ellip = self.lens.Lightlist[index]['ellipticity']
            ar = 1 - ellip

            lens_x_pix_list.append(lens_x_pix)
            lens_y_pix_list.append(lens_y_pix)
            mag_list.append(mag)
            re_list.append(re)
            n_list.append(n)
            ar_list.append(ar)
            pa_list.append(pa)

        return (lens_x_pix_list, lens_y_pix_list, mag_list, re_list,\
                n_list, ar_list, pa_list)

    def image_info(self, img_size, pix_scale):
        image_x, image_y, image_mu, image_dt = self.image_dat

        image_x_pix_list = (image_x + img_size/2)/pix_scale + 1
        image_y_pix_list = (image_y + img_size/2)/pix_scale + 1

        mag_list = []

        for index in range(self.source.n_light):
            mag = self.source.Lightlist[index]['Mag']
            mag_list.append(mag)

        mag_list = np.array(mag_list) -2.5 * np.log10(np.abs(image_mu))

        return (image_x_pix_list, image_y_pix_list, mag_list)

    def sim_image(self, img_size, pix_scale, psf_file, output_file,\
                 sky=0, workdir=dir_data+'/tmp', zeropoint=25.):
        image_x_pix, image_y_pix, image_mag =\
                self.image_info(img_size, pix_scale)
        lens_x_pix, lens_y_pix, lens_mag, lens_re, lens_n, lens_ar, lens_pa =\
                self.lens_info(img_size, pix_scale)

        img_size_pixel = int(img_size/pix_scale)

        cwd = os.getcwd()
        os.chdir(workdir)

        imgtmp = 'img.fits'
        fittmp = 'fit'
        psftmp = 'psf.fits'

        os.system('cp ' + psf_file + ' ' + psftmp)

        gf=galfit.galfit(imgtmp, fittmp, workdir+'/')
        gf.setpsf(psftmp)
        gf.setplatescale(pix_scale)
        gf.setzeropoint(zeropoint)
        gf.setboxsize(400)
        gf.psfovrsamp=1

        gf.setimageregion((0, img_size_pixel+1, 0, img_size_pixel+1))

        GalaxyComponentList = []
        for index1 in range(len(lens_x_pix)):
            GalaxyComponentList.append(\
                galfit.galfit_sersic(lens_x_pix[index1], lens_y_pix[index1],\
                    lens_mag[index1], Re=lens_re[index1], n=lens_n[index1],\
                    axis_ratio=lens_ar[index1], pa=lens_pa[index1],\
                    freelist=[False,False,False,False,False,False]))
            gf.addobject(GalaxyComponentList[-1])

        PSFComponentList = []
        for index2 in range(len(image_x_pix)):
            PSFComponentList.append(\
                galfit.galfit_psf(image_x_pix[index2], image_y_pix[index2],\
                    image_mag[index2], freelist=[False,False]))
            gf.addobject(PSFComponentList[-1])

        gf.addobject(galfit.galfit_sky(sky, freelist=[False]))
        gf.fit()

        os.chdir(cwd)
        os.system('mv ' + workdir + '/fit.fits ' + output_file)

def fit_galfit(psinfo, image, output_file, psf_file, pix_scale, size,
              workdir=dir_data+'/galfit/tmp', zeropoint=25.):
    # get the lensed images: pos & mag

    data = np.loadtxt(psinfo)

    image_x, image_y, image_mu, image_dt = data[1:].T
    print(image_x)

    # coordinate conversion
    x_size, y_size = size
    image_x_pix = image_x/pix_scale + 1 + x_size/2
    image_y_pix = image_y/pix_scale + 1 + y_size/2

    mag0 = 8.5
    image_mag = - 2.5 * np.log10(np.abs(image_mu))

    # sort it 
    sortidx = np.argsort(image_mag)
    image_x_pix = image_x_pix[sortidx]
    image_y_pix = image_y_pix[sortidx]
    image_mag = image_mag[sortidx]
    image_mag = image_mag - image_mag[0] + mag0

    cwd = os.getcwd()
    os.chdir(workdir)
    print(workdir)

    imgtmp = os.path.basename(image)
    fittmp = 'fit'
    psftmp = 'psf.fits'

    os.system('cp ' + psf_file + ' ' + psftmp)

    gf=galfit.galfit(imgtmp, fittmp, './')
    gf.setpsf(psftmp)
    gf.setplatescale(pix_scale)
    gf.setzeropoint(zeropoint)
    gf.setboxsize(50)
    gf.psfovrsamp=1

    gf.setimageregion((0, x_size+1, 0, y_size+1))

    # make the constraint file

    f = open(workdir + '/constraint.txt','w')
    f.write(\
'''1_2_3_4 x offset
1_2_3_4 y offset
1_2_3_4 mag offset
''')
    gf.setconstraints('constraint.txt')

#    GalaxyComponentList = []
#    for index1 in range(len(lens_x_pix)):
#        GalaxyComponentList.append(\
#            galfit.galfit_sersic(lens_x_pix[index1], lens_y_pix[index1],\
#                lens_mag[index1], Re=lens_re[index1], n=lens_n[index1],\
#                axis_ratio=lens_ar[index1], pa=lens_pa[index1],\
#                freelist=[False,False,False,False,False,False]))
#        gf.addobject(GalaxyComponentList[-1])

    PSFComponentList = []
    for index2 in range(len(image_x_pix)):
        if index2 == 0:
            PSFComponentList.append(\
            galfit.galfit_psf(image_x_pix[index2], image_y_pix[index2],\
                image_mag[index2], freelist=[True,True]))
        else:
            PSFComponentList.append(\
            galfit.galfit_psf(image_x_pix[index2], image_y_pix[index2],\
                image_mag[index2], freelist=[True, True]))

        gf.addobject(PSFComponentList[-1])

    gf.fit()

    os.chdir(cwd)
    os.system('mv ' + workdir + '/fit.fits ' + output_file)

def model_image_logL(x, image_file, psf_file, pixscale,
                    plotconfig=False):
    # parameters: f (scaling factor), xq, yq

    '''
    Step 1:
        Define the lensing system
    '''

    # fixed parameters
    # None

    # free params

    f, dxq, dyq, dxl, dyl, mag0 = x

    if np.abs(f-1)>1 or np.abs(dxq)>0.5 or np.abs(dyq)>0.5\
       or np.abs(dxl)>0.5 or np.abs(dyl)>0.5\
       or np.abs(mag0-8.5)>1:
        return -np.inf

    sigma_lens = 120 * np.sqrt(f)
    x_source = (dxq-0.076) * f
    y_source = (dyq+0.016) * f

    x_lens = 0
    y_lens = 0
    PA = 81.1527
    e = 1 - 0.6536

    # define the lens system
    lens_sie = {'Type':'sie', 'sigma':sigma_lens, 'r_core':0.0,\
               'ellipticity':e, 'pa':PA, 'x_center':x_lens, 'y_center':y_lens}

    lens_light = {'Type':'sersic', 'x_center':0., 'y_center':0., 'Mag':99,\
                  'ellipticity':0.3, 'pa':30, 're':0.7, 'n':1}
    source_light = {'Type':'point', 'x_center':x_source, 'y_center':y_source, 'Mag':20}

    lens = OneObject(0.5, mass_component=[lens_sie],\
               light_component=[])
    source = OneObject(5.1, light_component = [source_light])

    lenssystem = LensSystem(lens, source)
#    print(lenssystem.image_dat)

    if plotconfig:
        lenssystem.plot_lensconfig('./test.pdf',[-0.6,0.6], [-0.6,0.6])

    '''
    Step 2:
        coordinate conversion: glafic to galfit
    '''
    # center of mass of the NIFS image

    realimage = fits.open(image_file)[1].data
    com = center_of_mass(realimage)
    ycom, xcom = com

    # galaxy position in DECaLS
    xgal_decals = 0.3637084000000255 + dxl
    ygal_decals = 0.15028319999998488 + dyl

    # galaxy position in NIFS
    xgal_nifs = xgal_decals / pixscale + xcom + 1
    ygal_nifs = ygal_decals / pixscale + ycom + 1

    print(xgal_nifs, ygal_nifs)

    # image position in NIFS

    image_x, image_y, image_mu, image_dt = lenssystem.image_dat

    if len(image_x)<4:
        return -np.inf

    size = [60, 30]
    x_size, y_size = size
    image_x_pix = image_x / pixscale + xgal_nifs
    image_y_pix = image_y / pixscale + ygal_nifs

    image_mag = -2.5 * np.log10(np.abs(image_mu))
    image_mag = mag0 + image_mag - np.min(image_mag)

    '''
    Step 3:
        generate mock image
    '''

    # define files

    cwd = os.getcwd()

    workdir = dir_data + '/tmp/'
    os.chdir(workdir)

    imgtmp = 'image23_inverted.fits'
    fittmp = 'fit%d'%np.random.randint(10000000)
    psftmp = 'psf.fits'

#    os.system('cp ' + psf_file + ' ' + psftmp)
#    os.system('cp ' + image_file + ' ' + imgtmp)

    gf=galfit.galfit(imgtmp, fittmp, './')
    gf.setpsf(psftmp)
    gf.setplatescale(pixscale)
    gf.setzeropoint(25)
    gf.setboxsize(50)
    gf.psfovrsamp=1

    gf.setimageregion((0, 61, 0, 31))

    # make the constraint file

    f = open(workdir + '/constraint.txt','w')
    f.write(\
'''1_2_3_4 x offset
1_2_3_4 y offset
1_2_3_4 mag offset
''')
    gf.setconstraints('constraint.txt')

    PSFComponentList = []
    for index2 in range(len(image_x_pix)):
        PSFComponentList.append(\
            galfit.galfit_psf(image_x_pix[index2], image_y_pix[index2],\
                image_mag[index2], freelist=[False, False]))

        gf.addobject(PSFComponentList[-1])

    gf.fit()


    '''
    Step 4:
        get log likelihood
    '''

    mockheader = fits.open(workdir + '/' + fittmp + '.fits')[2].header
    # flux
#    factor = np.sum(realimage*mockimage) / np.sum(mockimage*mockimage)
#
#    residual = realimage - mockimage * factor

#    sigma = 13066.528
#    logL = - 0.5 * np.sum(residual**2 / sigma**2)

    chi2nu = mockheader['CHI2NU']
    logL = - 0.5 * chi2nu * 60 * 31

    print(logL)
    if not plotconfig:
        os.system('rm %s*'%fittmp)

    os.chdir(cwd)
    return logL

def run_mcmc():

    nwalkers = 20
    ndim = 6

    sampler = emcee.EnsembleSampler(nwalkers, ndim, model_image_logL,\
                                    args=(dir_data+'/tmp/image.fits',dir_data+'/tmp/psf.fits', 0.043))

    pos = [np.array([0.8,0,0,0,0, 8.8]) + 1e-2*np.random.randn(ndim) for i in range(nwalkers)]

    sampler.run_mcmc(pos, 4000, progress=True)

    pickle.dump(sampler.chain[:, 2000:, :].reshape((-1, ndim)), open( "save.p", "wb" ) )

def test_result():
    x_lens = 0
    y_lens = 0
    PA = 81.1527
    e = 1 - 0.6536

#    f = 0.9229179224994378
#    dxq = 0.01737233943102852
#    dyq = 0.004444596678861662

    f = 1
    dxq = 0
    dyq = 0

    sigma_lens = 120 * np.sqrt(f)
    x_source = (dxq-0.076) * f
    y_source = (dyq+0.016) * f

    print(x_source, y_source)

    lens_sie = {'Type':'sie', 'sigma':sigma_lens, 'r_core':0.0,\
               'ellipticity':e, 'pa':PA, 'x_center':0, 'y_center':0}

    lens_light = {'Type':'sersic', 'x_center':0., 'y_center':0., 'Mag':99,\
                  'ellipticity':0.3, 'pa':30, 're':0.7, 'n':1}
    source_light = {'Type':'point', 'x_center':x_source, 'y_center':y_source, 'Mag':20}

    lens = OneObject(0.5, mass_component=[lens_sie],\
               light_component=[])
    source = OneObject(5.1, light_component = [source_light])

    lenssystem = LensSystem(lens, source)
    lenssystem.plot_lensconfig('./test.pdf',[-0.6,0.6], [-0.6,0.6])

    psinfo = 'out_point.dat'
    output_file = './J0025_model.fits'
    raw_image = dir_data + '/J0025/galfit/image23.fits'
    image = dir_data + '/J0025/galfit/image23_inverted.fits'
    raw_psf_file = dir_data + '/J0025/galfit/Telluric.fits'
    psf_file = dir_data + '/J0025/galfit/Telluric_inverted.fits'
    pix_scale= 0.043

    workdir = dir_data + '/J0025/galfit'

    fit_galfit(psinfo, image, output_file, psf_file, pix_scale, size=[60,30],\
              workdir=workdir)


def main():
#    test_result()
#    image = fits.open(dir_data+'/J0025/galfit/fit.fits')[1].data

#    model_image_logL([1, 0, 0], dir_data+'/tmp/image.fits',\
#                    dir_data+'/tmp/psf.fits', 0.043)
#    run_mcmc()


    samples = pickle.load(open( "save.p", "rb" ))
    print(len(samples))
    fig = corner.corner(samples, labels=["$f$", "$dxq$", "$dyq$", "$dxl$", "$dyl$", "$mag0$"])
    fig.savefig("triangle.png")

    result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                             zip(*np.percentile(samples, [16, 50, 84],
                                                axis=0)))

#    print([e for e in result])
    values = [e[0] for e in result]

#    f = 0.9066748748251495
#    dxq = -0.013266853376475578
#    dyq = -0.0025425933537572805
#    dxl = 0.01977861918454816
#    dyl = -0.15900879744895358
#    mag0 = 8.494924644335699

    f, dxq, dyq, dxl, dyl, mag0 = values

    model_image_logL([f, dxq, dyq, dxl, dyl, mag0], dir_data+'/tmp/image.fits',\
                    dir_data+'/tmp/psf.fits', 0.043, plotconfig=True)


if __name__=='__main__':
    main()
