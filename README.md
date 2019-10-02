# Characterization-of-unfolding-of-a-protein
Using nonlinear curve fitting to fit data obtained from Circular Dichroism spectroscopy, and thereby studying the characteristics of proteins

CurveFit_Biophy.m - The function imports dataset Unf3.dat and permforms nonlinear regression by
making use of the lsqcurvefit() function.
It makes use of another function myExp.m which specifies the nonlinear
function, determined by temperature dependent thermodynamic parameters.
Finally it prints the estimate and the errors of the parameters

A detailed report is given in the PDF file
