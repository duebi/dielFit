# dielFit

MATLAB routine to estimate rates of gross production and community respiration from diel measurements of concentration of oxygen, particles, or other parameters. This routine is reported in a manuscript by Barone et al. (2019) and builds on the approach used by Nicholson et al. (2015). Estimates of community rates rely on the assumption of constant respiration throughout the day while the assumption on gross production varies based on three models:

1. Linear model, assuming constant production during daytime
2. Sinusoidal model, assuming that production scales linearly with light intensity
3. P vs. E model, including a parameterization for light saturation and photoinhibition

Rates estimates are obtained using linear least squares while rate uncertainties are obtained by bootstrapping the model residuals. Residual autocorrelation is tested using the Durbin-Watson test.

Requirements:
- Matlab ver 2013b or above
- Statistics and Machine Learning Toolbox
- The routine suncycle.m to compute solar elevation (by Ch. Begler, available at http://mooring.ucsd.edu/software/matlab/doc/toolbox/geo/suncycle.html)

Example of use in MATLAB’s Command Window:
x = 736497 + [0.1 0.25 0.4 0.7 1.1 1.25 1.35 1.65 1.8 2.05 2.15 2.45 2.6 2.9];
y = 200 + [0.121 -0.240 -0.223 0.259 -0.017 -0.051 -0.107 0.307 0.324 0.063 -0.026 -0.226 0.426 0.373];
[rates,residuals,varmat,fitline] = dielFit(x,y,22.75,-158,-10);

REFERENCES:

Barone, B., Nicholson, D. P., Ferrón, S., Firing, E., and Karl, D. M. (2019). The estimation of gross oxygen production and community respiration from autonomous time‐series measurements in the oligotrophic ocean. Limnology and Oceanography: Methods, 17(12), 650-664. doi:10.1002/lom3.10340

Nicholson, D. P., Wilson, S. T., Doney, S. C., and Karl, D. M. (2015). Quantifying subtropical North Pacific gyre mixed layer primary productivity from Seaglider observations of diel oxygen cycles. Geophysical Research Letters, 42(10), 4032-4039. doi:10.1002/2015GL063065.



Version r0.18 - September 2018

Benedetto Barone