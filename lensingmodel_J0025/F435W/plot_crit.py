import numpy as np
import matplotlib.pyplot as plt

dat = np.loadtxt('out_crit.dat')
print(dat.shape)

def plot_mu_thetaE():

    for i in range(dat.shape[0]):
    #    i0 = i
    #    i1 = i+1
    #    if i1==dat.shape[0]:
    #        i1=0

        crit_x = [dat[i][0], dat[i][4]]
        crit_y = [dat[i][1], dat[i][5]]
        caus_x = [dat[i][2], dat[i][6]]
        caus_y = [dat[i][3], dat[i][7]]

        plt.plot(crit_x, crit_y, 'r-')
        plt.plot(caus_x, caus_y, 'b-')

    plt.show()



def plot_lensconfig(number, x_range, y_range):


    image_dat = np.loadtxt('findimg_%d_point.dat'%number).T
    crit_dat = np.loadtxt('findimg_%d_crit.dat'%number).T

#    image_dat = self.image_dat
#    crit_dat = self.findcrit()

    image_x, image_y, image_mu, image_dt = image_dat
#    source_x = self.source.Lightlist[0]['x_center']
#    source_y = self.source.Lightlist[0]['y_center']

    with open('findsrc_%d_optresult.dat'%number, 'r') as f:
        while True:
            line = f.readline()
            if line.startswith('point  5.0700'):
                linebreaked = line.split(' ')
                source_x = float(linebreaked[4])
                source_y = float(linebreaked[5])
                break

#    print(linebreaked, source_x, source_y)
#    return 0

    xi1, yi1, xs1, ys1, xi2, yi2, xs2, ys2 = crit_dat
    print(crit_dat)

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
#    fig.savefig(prefix+'_config.pdf')
    plt.show()


if __name__=='__main__':
    plot_lensconfig(120, [-1,1], [-1,1])
