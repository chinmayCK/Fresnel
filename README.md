# Fresnel
This repository contains a useful tool for optics and photonics. Specifically, you will find matlab files to compute the Fresnel reflection and transmission coefficients for generic (multilayered) bianisotropic planar media. 
Please refer to the documentation.pdf in the repository for details regarding how to compute these coefficients. 


If you are using it for your research, you can cite the following published work. 
'Thermal spin photonics in the near-field of nonreciprocal media', C.Khandekar and Z.Jacob, New J. Phys. 21, 103030 (2019).
https://iopscience.iop.org/article/10.1088/1367-2630/ab494d/pdf
Thank you. Hope this tool is useful for your research or project work.  


We recently used this tool and discovered interesting relations between Fresnel coefficients that are applicable 
for many material classes. This led us to the discovery of new spin-resolved Kirchhoff's laws of thermal radiation. 
If you are interested in rediscovering these relations, follow the steps below:
1. Download 'verifySKL.zip' and unzip it. 
2. Go to the folder in matlab and run the command 'figure; publish('make_report.m','pdf')' in the command prompt.
   This will produce reflection and transmission coefficients for many materials described in the script 'make_report.m'
   and for user-specified parameters in the same. You can change these parameters later on. The generated file will look
   like 'make_report.pdf' uploaded in this repository. 
3. The above command generates multiple figures which are then published in the pdf. Please do not play around with the 
   figures (like select, save tools) while the code is generating it because it can cause unwanted arrangement in the 
   output pdf. On my computer with matlab 2018, this code takes only 9minutes to produce this report. 
