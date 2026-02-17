# J0025_magnification
Code for determining the lensing magnification of J0025-0145 (Panebianco+2026, DOI 10.3847/2041-8213/ae40aa). 

The files are organized as follows:

- galfitmodel_J0025: the image cutouts and Galfit .feedme files used to fit quasar and foreground galaxy data in the JWST NIRCam F070W and F480M filters
- lensingmodel_J0025: the files used to compute the maximum possible magnification for each image
- psfcompilation_J0025: the code used to interpolate the PSF models for the JWST NIRCam F070W and F480M filters at the location of the quasar
