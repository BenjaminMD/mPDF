"""
This module contains the mPDF function for the total PDF calculator.

The mPDF function calculates the "unnormalized" magnetic pair 
distribution function (mPDF) based on provided scale factors 
and correlation length.

Attributes:
    None

Todo:
    * dcalc is the PDF or just a scalar?
    * mc, mstruc needs to be refered to better use a class approach
        otherwise this will not work as intended or be slow
"""
from diffpy.mpdf import MagSpecies, MagStructure, MPDFcalculator
from diffpy.structure.parsers import getParser
import numpy as np


class MPDF_Wrapper:
    def __init__(self, phase_cif, config):
        self.phase_cif = phase_cif
        self.config = config
        self.rmax = config.R_val.rmax
        self.rmin = config.R_val.rmin
        self.qdamp = config.PDF.qdamp
        self.species = type('', (), {})()
        self.get_struc()
        
    def get_struc(self):
        pcif = getParser('cif')
        struc = pcif.parseFile(self.phase_cif)
        self.struc = struc

    def add_magnetic_species(self, name, ffparamkey, strucIdxs, basisvecs, kvecs):
        mspec = MagSpecies(
            struc=self.struc,
            rmaxAtoms=self.rmax,
            ffparamkey=ffparamkey
        )
        mspec.strucIdxs = [strucIdxs] # determined from previous inspection of the unit cell
        mspec.basisvecs = np.array(basisvecs) # any vector perpendicular to the c axis will be fine
        mspec.kvecs = np.array(kvecs)
        setattr(self.species, name, mspec)
        
    def set_up_MagCalc(self):
        self.mstruc = MagStructure(rmaxAtoms=self.rmax)
        species = [s for s in dir(self.species) if not s.startswith("__")]
        for spec in species:
            mspec = getattr(self.species, spec)
            self.mstruc.loadSpecies(mspec)
        self.mstruc.makeAll()
        mc = MPDFcalculator(
            magstruc=self.mstruc,
            qdamp=self.qdamp,
            rmax=self.rmax,
            rmin=self.rmin
        )
        self.mc = mc

    def register_mPDF_in_Structure(self, fit):
        fit.recipe.PDF.registerFunction(self.mpdf)
        fit.recipe.PDF.setEquation(fit.recipe.PDF.getEquation() + "+ mpdf(ordscale, parascale, xi)")

        fit.pgs['Tb3Ni'].phase.addObserver(fit.recipe.PDF.ordscale.notify)

        fit.recipe.addVar(fit.recipe.PDF.parascale, 0.0, tag='mag')
        fit.recipe.addVar(fit.recipe.PDF.ordscale, 0.0, tag='mag')
        fit.recipe.addVar(fit.recipe.PDF.xi, 0.0, tag='mag')
        
        fit.recipe.restrain(fit.recipe.PDF.ordscale,lb=0,ub=10.0,sig=0.0001)
        fit.recipe.restrain(fit.recipe.PDF.parascale,lb=0,ub=10.0,sig=0.0001) # set reasonable bounds
        fit.recipe.restrain(fit.recipe.PDF.xi,lb=0.5,ub=1000.0,sig=0.0001)
        

    def mpdf(self, ordscale, parascale, xi):
        """Calculate the "unnormalized" mPDF.

        Calculate the "unnormalized" magnetic pair distribution 
        function based on scale factors and correlation length.

        Args:
            ordscale (type): Ordered scale factor (mc.ordScale).
            parascale (type): Paramagnetic scale factor (mc.paraScale).
            xi (type): Correlation length (mstruc.corrLength).

        Returns:
            type: The "unnormalized" mPDF.

        Note:
            The function regenerates the spins in each call to use 
            the latest structural parameters.
        """
        self.mc.ordScale = ordscale
        self.mc.paraScale = parascale
        self.mstruc.corrLength = xi
        self.mstruc.makeAll() 
        rcalc, fcalc, dcalc = self.mc.calc(both=True)
        return dcalc
