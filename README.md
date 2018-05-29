# dielFit

MATLAB routine to estimate gross oxygen production (GOP) and Community Respiration (CR) from diel measurements of oxygen concentration in aquatic environments. 

Requirements:
- Matlab ver 2013b or above
- Statistics and Machine Learning Toolbox
- The routine suncycle.m to compute solar elevation (by Ch. Begler, available at http://mooring.ucsd.edu/software/matlab/doc/toolbox/geo/suncycle.html)

Both auxiliary/ and data/ folders needs to be added to your MATLAB path before calling mesoHot.m or other routines.
The auxiliary/ folder contains routines used within mesoHot.m
The data/ folder contains altimetric data and HOT hydrographic and biogeochemical observations.

Example of use (in MATLAB’s Command Window):

x = 736497 + [0.1 0.25 0.4 0.7 1.1 1.25 1.35 1.65 1.8 2.05 2.15 2.45 2.6 2.9];
y = 200 + [0.121 -0.240 -0.223 0.259 -0.017 -0.051 -0.107 0.307 0.324 0.063 -0.026 -0.226 0.426 0.373];
[rates,residuals] = dielFit(x,y,22.75,-158,-10);

Benedetto Barone - May 2018
