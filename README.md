# Fresnel
This repository contains a useful tool for optics and photonics. Specifically, you will find matlab files (functions) to compute the Fresnel reflection and transmission coefficients for generic (multilayered) bianisotropic planar media. 

You can refer to the documentation.pdf in the repository for details regarding how to compute these coefficients. 

If you are using it for your research, you can cite one of the following published works.

'Thermal spin photonics in the near-field of nonreciprocal media', C.Khandekar and Z.Jacob, New J. Phys. 21, 103030 (2019).
https://iopscience.iop.org/article/10.1088/1367-2630/ab494d/pdf

'New spin-resolved thermal radiation laws for nonreciprocal bianisotropic media', C.Khandekar, F.Khosravi, Z.Li and Z.Jacob, New J. Phys. 21, 123005 (2020). 
https://iopscience.iop.org/article/10.1088/1367-2630/abc988/pdf


If you are interested in validating the spin-resolved Kirchhoff's laws derived in the second work above, follow the steps below: 

1. Download the repository and go to the folder in matlab. 
2. Run the command 'figure; publish('make_report.m','pdf')' in the command prompt. 
   This will produce reflection and transmission coefficients for many materials described in the script 'make_report.m'
   and for user-specified parameters in the same. You can change these parameters later on. The generated file will look
   like 'make_report.pdf' uploaded in this repository. 
   
The above command generates multiple figures (around 20) each containing table of reflection and transmission coefficients. These figures are then published in the pdf. The code takes less than 1min to finish but capturing and publishing figures as a report can take some time (tpyically less than 3minutes). 

Thank you. Hope this tool is useful for your research or project work.  
