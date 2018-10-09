from astropy.io import ascii, fits
import glob
import os

fileList = glob.glob('new_k1255_ascii_w_halpha/w*.txt')

for oneFile in fileList:
    if oneFile != 'readme.txt':
        dat = ascii.read(oneFile)
        dat['col1'].name = 'Wavelength'
        dat['col2'].name = 'Flux'
        if len(dat.colnames) >= 3:
            dat['col3'].name = 'RawFlux'
        basename = os.path.splitext(os.path.basename(oneFile))[0]
        dat.write('new_k1255_fits_w_halpha/'+basename+'.fits')
    