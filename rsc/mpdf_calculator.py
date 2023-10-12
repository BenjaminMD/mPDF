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


def mpdf(mc, mstruc, ordscale, parascale, xi):
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
    mc.ordScale = ordscale
    mc.paraScale = parascale
    mstruc.corrLength = xi
    mstruc.makeAll() 
    rcalc, fcalc, dcalc = mc.calc(both=True)
    return dcalc
