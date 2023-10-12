## mPDF Module README

### Overview:
This module is an extension of `ezfit` and aims to incorporate the magnetic PDF fit into the fitting process. It contains the mPDF function for the total PDF calculator.

### Description:
The mPDF function calculates the "unnormalized" magnetic pair distribution function (mPDF) based on provided scale factors and correlation length.

### Attributes:
- None

### Dependencies:
- `ezfit`: For fitting PDF data.
- `diffpy.mpdf`: For magnetic species, magnetic structure, and MPDF calculator.
- `diffpy.structure.parsers`: For parsing CIF files.
- `numpy`: For numerical operations.
- `matplotlib`: For plotting.

### Class: MPDF_Wrapper
This class serves as a wrapper for the mPDF calculation. It initializes with a phase CIF and configuration, and provides methods to:
- Parse the CIF file and get the structure.
- Add magnetic species.
- Set up the magnetic calculator.
- Register the mPDF in the structure.
- Calculate the "unnormalized" mPDF.

### Usage:
1. Initialize the `FitPDF` with data and contributions.
2. Create an instance of `MPDF_Wrapper` with the CIF file and configuration.
3. Use the `add_magnetic_species` method to add magnetic species to the wrapper.
4. Set up the magnetic calculator using `set_up_MagCalc`.
5. Update the recipe using `fit.update_recipe()`.
6. Register the mPDF in the structure using `register_mPDF_in_Structure`.
7. Load results from a file using `fit.LoadResFromFile`.
8. Run the fit using `fit.run_fit()`.

### To-Do:
- understand add magentic species kvecs and basis and index
