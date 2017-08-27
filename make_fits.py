from astropy.io import ascii, fits
import glob
import os

fileList = glob.glob('kic1255_ascii/w*.txt')

for oneFile in fileList:
    if oneFile != 'readme.txt':
        dat = ascii.read(oneFile)
        dat['col1'].name = 'Wavelength'
        dat['col2'].name = 'Flux'
        if len(dat.colnames) >= 3:
            dat['col3'].name = 'RawFlux'
        basename = os.path.splitext(os.path.basename(oneFile))[0]
        dat.write('kic1255_fits/'+basename+'.fits')
    