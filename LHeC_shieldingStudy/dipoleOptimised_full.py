"""
C. Monaghan 11/08/2021

This script can be imported to allow access to various functions used in the
study of the LHeC IR regarding synchrotron radiation shielding.

## NOTICE ##
This lattice setup is for the SIMPLE IR (Full).
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

    def __init__(self, colMat, ngenerate, nruns, thickness, runKey):
        self._colMat    = colMat        # (the material of the shielding being studied)
        self._ngenerate = ngenerate     # (num primary particles i.e. eletrons)
        self._nruns     = nruns         # (num iterations)
        self._thickness = thickness     # metre
        self._buffer    = None          # (buffer containing the percentage from each run)
        self._runKey    = runKey

        self._seed     = 12             # a particular seed, for reproducability
        self.eAperture = 0.005          # metre (half x-y size)
        self.pAperture = 0.02           # metre (half x-y size)
        self.sep       = 0.121896       # metre (seperation of beam centroid)

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

    # Store cut data to compare the data which is lost at each stage
    def _storeZpCut(self, sampler, after=0, _partID=22, _xpos=None):
        self._cutDataBuffer.append([_np.asarray(len(sampler.data['energy'][(sampler.data['partID']==_partID)&(sampler.data['zp']>=0)])), len(sampler.data['energy'][sampler.data['partID']==22])])
    # end _storeCutData (func)
    
    def _getNum(self, sampler, _partID=22):
        """
        Return the number of photons (by default, can be any particle) at a particular sampler.

        By default what is returned is the number of photons at a sampler excluding photons in
        the electron or proton aperture.


        TODO:   Add functionality to turn of cuts and get number of photons in any regions, eg 
                for when you want to know the percentage of all photons passing through the aperture.
                This could be done with just a seperate function as it only needs to be called once,
                or maybe more if an average is to be calculated along with each study.
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
        a.AddIncludePre("extra-{}.gmad".format(self._runKey))

        # Start definition of lattice
        a.AddDrift('DRIFT_0', 5)
        a.AddDipole('BEND_0', length=20, angle=0.0244)
        #a.AddDipole('BEND_1', length=10, angle=0.0122)
        a.AddDrift('DRIFT_1', 5-self._thickness)

        # Definitions and placements of the shielding material being studied.
        # Do not reccomend changing this unless required, the material placement is determined 
        # based in the definitions given at the start.
        width = (self.sep-self.eAperture)*2 # width of first collimater so that proton aperture is in correct place
        a.AddRCol('COL_0', self._thickness, material=self._colMat, xsize=self.pAperture, ysize=self.pAperture, horizontalWidth=width, offsetX=((width/2)+self.eAperture))
        a.AddPlacement('COL_1_p', bdsimElement='"COL_1"', referenceElement='"COL_0"', x=((width/2)+width+self.eAperture))
        # open the extra.gmad file and replace it with new definition with identical thickness and material
        f = open("GMAD/extra-{}.gmad".format(self._runKey), "w")
        f.write('COL_1: rcol, horizontalWidth={}, l={}, material="{}", xsize=0.0, ysize=0.0;'.format(width,self._thickness, self._colMat))
        f.close()

        # add samplers at the end of each element in the lattice
        a.AddSampler('all')

        # Begin beam definition
        b = _pybdsim.Beam.Beam(X0=0.0,
        	                Xp0=0.0,
        	                Y0=0.0,
        	                Yp0=0.0,
        	                alfx=150,
        	                alfy=150,
        	                betx=2250,
        	                bety=2250,
        	                dispx=0.1337,
        	                dispxp=0.0121997,
        	                dispy=0.1337,
        	                dispyp=0.0121997,
        	                distrtype="gausstwiss",
        	                emitx=0.5e-09,
        	                emity=0.5e-09,
        	                energy=50,
        	                particle="e-",
        	                sigmaE=0.00028)
        a.AddBeam(b)

        # Begin definition of simulation options 
        o = _pybdsim.Options.Options(
        	                        magnetGeometryType='"none"',
        	                        preprocessGDML=0,
        	                        physicsList='"synch_rad em"',
        	                        hStyle=1, 
                                    beampipeMaterial='"Cu"',
                                    apertureType='"elliptical"', 
	                                aper1=0.5,
	                                aper2=0.3,
                                    horizontalWidth=1.05,
                                    worldMaterial='"vacuum"',
                                    maximumStepLength=0.1)
        a.AddOptions(o)

        # Path is relative to where run from so be careful these directories are created before 
        # the start of running, for example I have used the os package to ensure each run is in same place
        a.Write("GMAD/input-{}".format(self._runKey)) # Write the gmad output to this location.
    # end genGMAD (func)

    def genRebdsim(self):
        """
        Generate a rebdsim file containing all the histograms which are required to produce the desired output data.

        See BDSIM docs for rebdsim examples and explanations. 
        """
        f = open("rebdsim-input.txt")
        lines = ["SimpleHistogram1D Event. NPhotons_DRIFT_1_cuts_1 {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22&DRIFT_1.zp>=0&DRIFT_1.x>={}&DRIFT_1.x<={}".format(self.eAperture, (self.sep-self.pAperture)), 
                "SimpleHistogram1D Event. NPhotons_DRIFT_1_cuts_2 {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22&DRIFT_1.zp>=0&DRIFT_1.x>={}".format(self.sep+self.pAperture),
                "SimpleHistogram1D Event. NPhotons_COL_0_cuts_1 {{100}} {{0:0.8}} COL_0.x COL_0.partID==22&COL_0.zp>=0&COL_0.x>={}&COL_0.x<={}".format(self.eAperture, (self.sep-self.pAperture)), 
                "SimpleHistogram1D Event. NPhotons_COL_0_cuts_2 {{100}} {{0:0.8}} COL_0.x COL_0.partID==22&COL_0.zp>=0&COL_0.x>={}".format(self.sep+self.pAperture),
                "SimpleHistogram1D Event. NPhotons_eAper {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22&DRIFT_1.zp>=0&DRIFT_1.x<={}&DRIFT_1.x>={}".format(self.eAperture, (-self.eAperture)),
                "SimpleHistogram1D Event. NPhotons_pAper {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22&DRIFT_1.zp>=0&DRIFT_1.x>={}&DRIFT_1.x>={}".format(self.sep-self.pAperture, (self.sep+self.pAperture)),
                "SimpleHistogram1D Event. NPhotons_DRIFT_1_total {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22",
                "SimpleHistogram1D Event. NPhotons_DRIFT_1_zp {{100}} {{0:0.8}} DRIFT_1.x DRIFT_1.partID==22&DRIFT_1.zp>=0"]
        f.writelines(lines)

    def runStudy(self):
        """
        Runs the set study, must call genGMAD() before running (unless provided files manually).
        If providing manually the main gmad must be in the directory as 'GMAD/input.gmad'.

        The outputs provide the average over a set of nruns each with a different seed.
        By default the before and after samplers are named DRIFT_0 and COL_0 respectively.
        
        Also returned is the standard error on the value. + the range as an array with two values
        """
        #self.genRebdsim()
        _buffer = [] # buffer to store the percentage absorbed on each run 
        for i in range(self._nruns):
            # Imprtant to make sure this directory exists where python being called from
            outfile = 'DATA/{}/{}-{}m'.format(self._runKey,self._colMat,self._thickness)

            runOptions = "--seed={}".format((i*42)+23)

            # run bdsim
            _pybdsim.Run.Bdsim("GMAD/input-{}.gmad".format(self._runKey), outfile, ngenerate=self._ngenerate, options=runOptions)

            _pybdsim.Run.Rebdsim("rebdsim-input.txt", "{}.root".format(outfile), "tmp/rebdsim-{}.root".format(self._runKey))
            
            # load the bdsim data from last run
            d = _pybdsim.Data.Load("tmp/rebdsim-{}.root".format(self._runKey))

            self._totalPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_DRIFT_1_total'].entries)
            self._eAperPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_eAper'].entries)
            self._pAperPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_pAper'].entries)
            self._zpPhotons.append(d.histogramspy['Event/SimpleHistograms/NPhotons_DRIFT_1_zp'].entries)

            numBefore = (d.histogramspy['Event/SimpleHistograms/NPhotons_DRIFT_1_cuts_1'].entries + d.histogramspy['Event/SimpleHistograms/NPhotons_DRIFT_1_cuts_2'].entries)
            numAfter  = (d.histogramspy['Event/SimpleHistograms/NPhotons_COL_0_cuts_1'].entries + d.histogramspy['Event/SimpleHistograms/NPhotons_COL_0_cuts_2'].entries)
            
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