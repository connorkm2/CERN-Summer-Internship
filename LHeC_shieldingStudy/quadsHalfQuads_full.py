"""
C. Monaghan 11/08/2021

This script can be imported to allow access to various functions used in the
study of the LHeC IR regarding synchrotron radiation shielding.

## NOTICE ##
This lattice setup is for the optimised IR (Full).
Change seperation when other IR is used (self.sep).

TODO:   Add a option to write data at end of study, can prevent complete data
        loss when a crash occurs.
"""

import numpy as _np
from numpy.core.fromnumeric import std
import pybdsim as _pybdsim
import sys

class shieldingStudy:
    """
    Run a study of a particular material for a particular thickness.

    Example:

    >>> s = shieldingStudy(colMat="Cu", ngenerate=10000, nruns=10, thickness=0.05)
    >>> s.genGMAD()
    >>> value, err = s.runStudy()

    Optional parameters can be set afterwards like changing the electron aperture
    and the proton aperture including the seperation of the two beam centroids.  
    """

    def __init__(self, colMat, ngenerate, nruns, thickness, runKey, extraShielding=False, extraT=0):
        self._colMat    = colMat        # (the material of the shielding being studied)
        self._ngenerate = ngenerate     # (num primary particles i.e. eletrons)
        self._nruns     = nruns         # (num iterations)
        self._thickness = thickness     # metre
        self._buffer    = None          # (buffer containing the percentage from each run)
        self._runKey    = runKey

        self._seed     = 12             # a particular seed, for reproducability
        self.eAperture = 0.00015        # metre (half x-y size)
        self.pAperture = 0.02           # metre (half x-y size)
        self.sep       = 0.106105       # metre (seperation of beam centroid)
        self._colNames = ["COL_END", "COL_BEND"]    # (names of collimaters)
        self._extra    = extraShielding
        self._extraT   = extraT

        # Storage for particular data
        self._totalPhotons = []
        self._eAperPhotons = []
        self._pAperPhotons = []
        self._zpPhotons    = []
        
        ## Unfinished: store cut data to compare the data which is lost at each stage
        #self._cutDataBuffer = []
        #self._cutCounter    = 0
        self._cutDataBuffer = []
        self._aperBuffer    = []
    # end __init__ (func)

    def _storeZpCut(self, sampler, after=0, _partID=22, _xpos=None):
        """
        Store cut data to compare the data which is lost at each stage
        """
        self._cutDataBuffer.append([_np.asarray(len(sampler.data['energy'][(sampler.data['partID']==_partID)&(sampler.data['zp']>=0)])), len(sampler.data['energy'][sampler.data['partID']==22])])
    # end _storeCutData (func)
    
    def _getNum(self, sampler, _partID=22):
        """
        Return the number of photons (by default, can be any particle) at a particular sampler.

        By default what is returned is the number of photons at a sampler excluding photons in
        the electron or proton aperture.
        """
        # Calc number of particles with _partID impinging between two apertures
        num_1   = len(sampler.data['energy'][(sampler.data['partID']==_partID) 
                                            &(sampler.data['zp']>=0)
                                            &(sampler.data['x']>=self.eAperture)
                                            &(sampler.data['x']<=(self.sep-self.pAperture))])
        # Calc the rest of particles with _partID impinging beyond p aperture
        num_2   = len(sampler.data['energy'][(sampler.data['partID']==_partID)
                                            &(sampler.data['zp']>=0)
                                            &(sampler.data['x']>=(self.sep+self.pAperture))])
        return num_1+num_2 
    # end _getNum (func)

    def _getNumAper(self, sampler1, _partID=22, sampler2=None, studyAfter=False):
        """
        Get the number of photons passing through the apertures. 

        To study the number of photons after the aperture turn on 'studyAfter'.
        This may be useful if one would like to investigte the number of photons 
        interacting inside the aperture of the collimator.

        'sampler2' is therefore the same type as 'sampler1' but for after the aperture. 
        """
        # Calc number of particles with _partID passing through e aperture
        eAper   = len(sampler1.data['energy'][(sampler1.data['partID']==_partID) 
                                            &(sampler1.data['zp']>=0)
                                            &(sampler1.data['x']<=self.eAperture)
                                            &(sampler1.data['x']>=(-self.pAperture))])
        # Calc the rest of particles with _partID passing through p aperture
        pAper   = len(sampler1.data['energy'][(sampler1.data['partID']==_partID)
                                            &(sampler1.data['zp']>=0)
                                            &(sampler1.data['x']>=(self.sep-self.pAperture))
                                            &(sampler1.data['x']<=(self.sep+self.pAperture))])
        if studyAfter: 
            # Calc number of particles with _partID passing through e aperture
            eAper_a = len(sampler2.data['energy'][(sampler2.data['partID']==_partID) 
                                                &(sampler2.data['zp']>=0)
                                                &(sampler2.data['x']<=self.eAperture)
                                                &(sampler2.data['x']>=(-self.pAperture))])
            # Calc the rest of particles with _partID passing through p aperture
            pAper_a = len(sampler2.data['energy'][(sampler2.data['partID']==_partID)
                                                &(sampler2.data['zp']>=0)
                                                &(sampler2.data['x']>=(self.sep-self.pAperture))
                                                &(sampler2.data['x']<=(self.sep+self.pAperture))])
            return [eAper, pAper, eAper_a, pAper_a]
        else: 
            return [eAper, pAper]

    def _addExtraShielding(self, a):
        sep = 0.029
        eAp = 0.0025
        pAp = 0.01035/2
        width = (sep-eAp)*2 # width of first collimater so that proton aperture is in correct place
        a.AddRCol('{}_0'.format(self._colNames[1]), self._extraT, material=self._colMat, xsize=pAp, ysize=pAp, horizontalWidth=width, offsetX=((width/2)+eAp))

    def genGMAD(self):
        """
        Generate a set of GMAD files to the particular specification of this study as defined
        by the passed parameters when intiating this instance.

        This function is important to change when studying a different IR. 

        The shielding material parameters are calculated based on the aperture sizes.

        TODO:   This could be adapted to allow the user to pass there own machine/options/beam
                resulting in no need to adapt this function.
                this could be more streamlined and allow for much more studies
        """

        a = _pybdsim.Builder.Machine()

        # Add two gmad files which contain extra information. The first contains different 
        # definitions of concrete (reccomend to comment out this line if not required). 
        # The second file is essential for the collimation to work. The extra.gmad file contains 
        # the definition of extra externally placed collimater which ensures the whole synchrotron 
        # fan is incident on the material.
        a.AddIncludePre("material_Concretes.gmad")
        a.AddIncludePre("extra.gmad")

        # Start definition of lattice
        a.AddDrift('uDRIFT50', 0.5)
        a.AddDrift('uDRIFT_Q0', 1.871978)
        a.AddDrift('uDRIFT30_0', 0.3)
        a.AddDipole('uBEND_QY', angle=0.002072436, length=2.1771258, k1=-18.36839/166.7778)   
        a.AddDrift('uDRIFT20', 0.2)
        a.AddDipole('uBEND_QX', angle=0.002072436, length=2.1771258, k1=30.13370/166.7778)    
        a.AddDrift('uDRIFT30_1', 0.3-self._extraT)
        if self._extra: self._addExtraShielding(a)  # extra shielding to reduce synchrotron fan
        a.AddDipole('BEND_0', length=7.85153*2, angle=0.007473978*2)    # long dipole centre is IP
        a.AddDrift('dDRIFT30_0', 0.3)
        a.AddDipole('dBEND_QX', angle=0.002072436, length=2.1771258, k1=30.13370/166.7778)    
        a.AddDrift('dDRIFT20', 0.2)
        a.AddDipole('dBEND_QY', angle=0.002072436, length=2.1771258, k1=-18.36839/166.7778)  
        a.AddDrift('dDRIFT30_1', 0.3)
        a.AddDrift('dDRIFT_Q0', 1.871978)
        a.AddDrift('dDRIFT50', 0.5-self._thickness)

        ##
        # Definitions and placements of the shielding material being studied.
        # Do not reccomend changing this unless required, the material placement is determined 
        # based in the definitions given at the start.
        width = (self.sep-self.eAperture)*2 # width of first collimater so that proton aperture is in correct place
        a.AddRCol('{}_0'.format(self._colNames[0]), self._thickness, material=self._colMat, xsize=self.pAperture, ysize=self.pAperture, horizontalWidth=width, offsetX=((width/2)+self.eAperture))
        a.AddPlacement('{}_1_p'.format(self._colNames[0]), bdsimElement='"{}_1"'.format(self._colNames[0]), referenceElement='"{}_0"'.format(self._colNames[0]), x=((width/2)+width+self.eAperture))
        a.AddPlacement('{}_2_p'.format(self._colNames[0]), bdsimElement='"{}_2"'.format(self._colNames[0]), referenceElement='"{}_0"'.format(self._colNames[0]), x=((width/2)+(width*2)+self.eAperture))
        # open the extra.gmad file and replace it with new definition with identical thickness and material
        f = open("GMAD/extra.gmad", "w")
        f.write('{}_1: rcol, horizontalWidth={}, l={}, material="{}", xsize=0.0, ysize=0.0;\n {}_2: rcol, horizontalWidth={}, l={}, material="{}", xsize=0.0, ysize=0.0;\n {}_1: rcol, horizontalWidth={}, l={}, material="{}", xsize=0.0, ysize=0.0;'.format(self._colNames[0],width,self._thickness, self._colMat,self._colNames[0],width,self._thickness, self._colMat, self._colNames[1],(0.029-0.005)*2,self._extraT, self._colMat))
        f.close()
        ##

        # add samplers at the end of each element in the lattice
        a.AddSampler('all')

        # Begin beam definition
        b = _pybdsim.Beam.Beam(X0=0.0,
	                           Xp0=0.0, 
	                           Y0=0.0, 
	                           Yp0=0.0, 
	                           alfx=-0.035622, 
	                           alfy=99.948526, 
	                           betx=0.090510, 
	                           bety=4881.193917, 
	                           dispx=0.1337, 
	                           dispxp=0.0121997, 
	                           dispy=0.1337, 
	                           dispyp=0.0121997, 
	                           distrtype="gausstwiss", 
	                           emitx=5e-10, 
	                           emity=5e-10, 
	                           energy=50, 
	                           particle="e-", 
	                           sigmaE=0.00028)
        a.AddBeam(b)

        # Begin definition of simulation options 
        o = _pybdsim.Options.Options(
        	                        magnetGeometryType='"none"',    # use double here so when written to gmad file it outputs with "..."
        	                        preprocessGDML=0,
        	                        physicsList='"synch_rad em"',
        	                        hStyle=1, 
                                    beampipeMaterial='"Cu"',
                                    apertureType='"elliptical"', 
	                                aper1=0.5,
	                                aper2=0.3,
                                    horizontalWidth=1.05,
                                    worldMaterial='"vacuum"',
                                    maximumStepLength=0.1,
                                    integratorSet='"geant4"')
        a.AddOptions(o)

        # Path is relative to where run from so be careful these directories are created before 
        # the start of running, for example I have used the os package to ensure each run is in same place
        a.Write('GMAD/input') # Write the gmad output to this location.
    # end genGMAD (func)

    def genRebdsim(self):
        """
        Generate a rebdsim file containing all the histograms which are required to produce the desired output data.

        See BDSIM docs for rebdsim examples and explanations. 
        """
        f = open("rebdsim-input.txt", "w")
        lines = ["SimpleHistogram1D Event. NPhotons_dDRIFT50_cuts_1 {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22&dDRIFT50.zp>=0&dDRIFT50.x>={}&dDRIFT50.x<={}".format(self.eAperture, (self.sep-self.pAperture)), 
                "\nSimpleHistogram1D Event. NPhotons_dDRIFT50_cuts_2 {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22&dDRIFT50.zp>=0&dDRIFT50.x>={}".format(self.sep+self.pAperture),
                "\nSimpleHistogram1D Event. NPhotons_COL_END_0_cuts_1 {{100}} {{0:0.8}} COL_END_0.x COL_END_0.partID==22&COL_END_0.zp>=0&COL_END_0.x>={}&COL_END_0.x<={}".format(self.eAperture, (self.sep-self.pAperture)), 
                "\nSimpleHistogram1D Event. NPhotons_COL_END_0_cuts_2 {{100}} {{0:0.8}} COL_END_0.x COL_END_0.partID==22&COL_END_0.zp>=0&COL_END_0.x>={}".format(self.sep+self.pAperture),
                "\nSimpleHistogram1D Event. NPhotons_eAper {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22&dDRIFT50.zp>=0&dDRIFT50.x<={}&dDRIFT50.x>={}".format(self.eAperture, (-self.eAperture)),
                "\nSimpleHistogram1D Event. NPhotons_pAper {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22&dDRIFT50.zp>=0&dDRIFT50.x>={}&dDRIFT50.x>={}".format(self.sep-self.pAperture, (self.sep+self.pAperture)),
                "\nSimpleHistogram1D Event. NPhotons_dDRIFT50_total {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22",
                "\nSimpleHistogram1D Event. NPhotons_dDRIFT50_zp {{100}} {{0:0.8}} dDRIFT50.x dDRIFT50.partID==22&dDRIFT50.zp>=0"]
        f.writelines(lines)
        f.close()

    def runStudy(self):
        """
        Runs the set study, must call genGMAD() before running (unless provided files manually).
        If providing manually the main gmad must be in the directory as 'GMAD/input.gmad'.

        The outputs provide the average over a set of nruns each with a different seed.
        By default the before and after samplers are named DRIFT_0 and self._colNames respectively.
        
        Also returned is the standard error on the value. + the range as an array with two values

        TODO: Needs updating to study extra material
        """
        self.genRebdsim()
        _buffer = [] # buffer to store the percentage absorbed on each run 
        for i in range(self._nruns):
            # Imprtant to make sure this directory exists where python being called from
            outfile = 'DATA/run_{}_{}/{}-{}m_{}'.format(self._colMat,self._runKey,self._colMat,self._thickness, self._runKey)

            runOptions = "--seed={}".format((i*42)+23)

            # run bdsim
            _pybdsim.Run.Bdsim('GMAD/input.gmad', outfile, ngenerate=self._ngenerate, options=runOptions)

            # run rebdsim
            _pybdsim.Run.Rebdsim("rebdsim-input.txt", 'DATA/run_{}_{}/{}-{}m_{}.root'.format(self._colMat,self._runKey,self._colMat,self._thickness, self._runKey), "tmp/rebdsim-{}.root".format(self._runKey))
  
            # load the bdsim data from last run
            d = _pybdsim.Data.Load("tmp/rebdsim-{}.root".format(self._runKey))

            self._totalPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_dDRIFT50_total'].entries)
            self._eAperPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_eAper'].entries)
            self._pAperPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_pAper'].entries)
            self._zpPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_dDRIFT50_zp'].entries)

            numBefore = (d.histogramspy['Event/SimpleHistograms/NPhotons_dDRIFT50_cuts_1'].entries + d.histogramspy['Event/SimpleHistograms/NPhotons_dDRIFT50_cuts_2'].entries)
            numAfter  = (d.histogramspy['Event/SimpleHistograms/NPhotons_COL_END_0_cuts_1'].entries + d.histogramspy['Event/SimpleHistograms/NPhotons_COL_END_0_cuts_2'].entries)
            
            # append to buffer the percentae absorbed
            _buffer.append(1-(numAfter/numBefore))
            self._buffer = _buffer
        
        # calculate the mean and the standard error
        value     = _np.mean(_np.asarray(_buffer))
        err       = (_np.std(_np.asarray(_buffer)))/(_np.sqrt(len(_buffer)))
        val_range = _np.asarray([min(_buffer), max(_buffer)])
        return value, err, val_range
    # end runStudy (func)

    def getBuffer(self):
        return self._buffer
    # end getBuffer (func)

    def getTotalPhotons(self):
        t = _np.asarray(self._totalPhotons)
        return _np.mean(t), _np.std(t)/_np.sqrt(self._nruns)
    # end getBuffer (func)

    def getTotalAper(self):
        e = _np.asarray(self._eAperPhotons)
        p = _np.asarray(self._pAperPhotons)
        return _np.mean(e), _np.mean(p), _np.std(e)/_np.sqrt(self._nruns), _np.std(p)/_np.sqrt(self._nruns)
    # end getBuffer (func)

    def getZpCut(self):
        z = _np.asarray(self._zpPhotons)
        return _np.mean(z), _np.std(z)/_np.sqrt(self._nruns)
    
# end shieldingStudy (class)