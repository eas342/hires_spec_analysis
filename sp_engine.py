from astropy.io import ascii, fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import es_gen
import pdb
import glob
from scipy.interpolate import interp1d
from astropy.convolution import Gaussian1DKernel, convolve
import os
import es_gen

class stellModel():
    def __init__(self,fileName,waveOffset=-2.0):
        """ A class to hold BOSZ stellar models """
        self.fileName = fileName
        HDUList = fits.open(fileName)
        self.header = HDUList[0].header
        self.logg = self.header['LOGG']
        self.temp = self.header['T_EFF']
        self.data = HDUList[1].data
        
        HDUList.close()
        self.Res = self.header['INSBROAD']
        self.waveOffset = waveOffset
        self.data['Wavelength'] = self.data['Wavelength'] + self.waveOffset
    
    def getPoints(self,wavRange):
        pts = (self.data['Wavelength'] > wavRange[0]) & (self.data['Wavelength'] < wavRange[1])
        return pts
    
    def getSpec(self,normRange=[5176,5180],showRange=[5160,5190],doSmooth=False,
                removePoly=None):
        """ Gets a spectrum and does normalization and removing of polynomials,
        if necessary
        normRange: 2 element list with wavelength range in AA
            The wavelength range over which to perform normalization
        showRange: 2 element list with wavelength range in AA
            The wavelength range over which to restrict the output array
        doSmooth: bool
            Use the smoothed spectrum self.yconv? (created with self.smooth())
        removePoly: bool
            Remove a polynomial from the spectrum self.yModel (created
            with self.removePoly?)
        """
        pts = self.getPoints(normRange)
        
        showPts = self.getPoints(showRange)
        if doSmooth == True:
            showY = self.yconv
        else:
            showY = self.data['SpecificIntensity']
        
        if removePoly == True:
            showY = showY / self.yModel
        
        norm = np.nanmedian(showY[pts])
        return self.data['Wavelength'][showPts], showY[showPts]/norm
        
    
    def plot(self,ax=None,offset=0.,customName=None,
             normRange=[5176,5180],showRange=[5160,5190],doSmooth=False,removePoly=None):
        
        x, y = self.getSpec(normRange=normRange,showRange=showRange,doSmooth=doSmooth,
                           removePoly=removePoly)
        if customName == None:
            doLabel = 'logg = '+str(self.logg)+',T='+str(self.temp)+' K'
        else:
            doLabel = customName
        
        ax.plot(x,y + offset,label=doLabel)
        
    def smooth(self,refWavel,outRes):
        """ Smooths at a reference wavelength to given resolution """
        self.outRes = outRes
        if outRes >= self.Res:
            raise ValueError('Output Res must be lower than intrinsic!')
            
        ## Calculate the convolution kernel
        self.RConv = 1./np.sqrt((1./outRes)**2 - (1./self.Res)**2)
        ind = np.searchsorted(self.data['Wavelength'],refWavel)
        dX = self.data['Wavelength'][ind+1] - self.data['Wavelength'][ind]
        stDevAngstroms = (refWavel/self.RConv) ## sigma in Angstroms
        stdDevConv = stDevAngstroms / dX ## in pixels
        self.modelKernel = Gaussian1DKernel(stddev=stdDevConv)
        self.yconv = convolve(self.data['SpecificIntensity'], 
                              self.modelKernel, boundary='extend')

    def removePoly(self,waveRange=[5176,5180],nPoly=1,**kwargs):
        """ Remove a polynomial trend to get a better normalized spec """
        pts = self.getPoints(waveRange)
        polyFit = es_gen.robust_poly(self.data['Wavelength'][pts],
                                    self.data['SpecificIntensity'][pts],nPoly,
                                    **kwargs)
        self.yModel = np.polyval(polyFit,self.data['Wavelength'])
        self.yDetrend = self.data['SpecificIntensity'] / self.yModel
        

class measuredArray(stellModel):
    def __init__(self,wavel,intensity,Res=1e5):
        """ A class for measured data that has the same methods as stellModel """
        self.fileName = 'from_array'
        self.logg= 'unknown'
        self.Res = Res
        self.data = Table()
        self.data['Wavelength'] = wavel
        self.data['SpecificIntensity'] = intensity
        self.yconv = None
        
    def chiSquared(self,model,waveRange=[5160,5190],removePoly=True,doSmooth=True,
                   normRange=[5176,5179],deltaWInterp = 1.,SNR=100.):
        """ Returns the chi-squared value over an interval
        Parameters:
        --------------------
        waveRange: 2 element list
            Wavelength region over which to calculate chi-squared
        deltaWInterp: float
            angstroms of extra tails to interpolate model onto spectrum
        Returns:
        --------------------
        chiSQ: float
            The chi-squared value
        chiString: str
            A string that provides a useful chi^2 for placing in a plot legend
        """
        modX, modY = model.getSpec(doSmooth=doSmooth,removePoly=removePoly,normRange=normRange,
                                   showRange=[waveRange[0]- deltaWInterp,waveRange[1] + deltaWInterp])
        fInterp = interp1d(modX,modY)
        
        dataX, dataY = self.getSpec(doSmooth=doSmooth,removePoly=removePoly,normRange=normRange,
                                    showRange=waveRange)
        
        modInterp = fInterp(dataX)
        roughChiSQ = np.sum((dataY - modInterp)**2 * SNR**2)
        
        chiString = r' $\chi^2$='+"{:8.0f}".format(roughChiSQ)
        return roughChiSQ, chiString
        
class subaruSpec():
    def __init__(self,fileName,waveOffset=0.0):
        """ Get a subaru spectrum and normalize 
        Parameters
        offset: float
            Put a wavelength offset - not sure why they're off
        """
        self.fileName = fileName
        self.baseName = os.path.basename(fileName)
        self.dat = Table(fits.getdata(fileName))
        self.dat.sort('Wavelength')
        fileNumText = os.path.splitext(self.baseName)[0].split('w')[-1]
        self.fileNumber = int(fileNumText)
        if np.mod(self.fileNumber,2) == 0:
            self.ccdSide = 'blue'
        else:
            self.ccdSide = 'red'
        
        self.waveOffset = waveOffset
        self.wave = self.dat['Wavelength'] + self.waveOffset
        self.normalize(region=[5160,5195])
        
    def normalize(self,region=[5160,5195]):
        closepts = (self.wave > 5160.) & (self.wave < 5195.)
        normValue1 = np.median(self.dat['Flux'][closepts])
        self.normFlux = self.dat['Flux'] / normValue1
        
    def plot(self,ax=None,offset=0.0):
        ax.plot(self.wave,self.normFlux + offset,label=self.baseName)


