function [rates,residuals,varmat] = dielFit(mtime,conc,lat,lon,varargin)%EmaxEk,EmaxEb,alpha)
%
% function [rates,residuals,varmat,fitline] = dielFit(mtime,conc,lat,lon,timezone)
%          [rates,residuals,varmat,fitline] = dielFit(mtime,conc,lat,lon,timezone,EmaxEk,EmaxEb)
%          [rates,residuals,varmat,fitline] = dielFit(mtime,conc,lat,lon,timezone,EmaxEk,EmaxEb,alpha)
%
%   Computes rates of gross production and community respiration from diel
%   cycles of oxygen concentration thorough ordinary linear least squares.
%   Parameter confidence intervals and the variance-covariance matrix are
%   computed via bootstrapping of residual with 200 iterations
%   Residual autocorrelation is tested using the Durbin-Watson test.
%
%   Three models are used for rate estimation:
%       1. Linear model - constant production during daytime
%       2. Sinusoidal model - production scales linearly with light intensity
%       3. P vs. E model - includes light saturation and photoinhibition
%
%   This function requires MATLAB ver 2013b or above, the Statistics and
%   Machine Learning Toolbox, and the routine suncycle.m to compute solar elevation
%   (by C. Begler, available at http://mooring.ucsd.edu/software/matlab/doc/toolbox/geo/suncycle.html)
%
% INPUT:
%   dates: local date and time in Matlab format (datenum)
%   conc: concentration of oxygen or particles at each value of dates
%   lat: latitude in degrees North (Station ALOHA = 22.75)
%   lon: longitude in degrees East (Station ALOHA = -158)
%   timezone: hours from GMT time to local time (Station ALOHA = -10)
%   EmaxEk (optional): light saturation parameter, maximum daily irradiance divided by Ek (DEFAULT = 1)
%   EmaxEb (optional): photoinhibition parameter (DEFAULT = 0)
%   alpha (optional): defines level of percent confidence intervals as 100*(alpha-1) (DEFAULT = 0.317)
% OUTPUT:
%   rates: a table containing the values of gross production and
%       respiration for the three models plus fit statistics (coefficient
%       of determination, p-value, (1-alpha)x100% confidence intervals for
%       GPP and R, p-value from the Durbin-Watson test).
%       Rate units are the same as input variable 'conc' per day.
%   residuals: a matrix with the residuals from the three models reported
%       in separate rows: linear on the 1st row, sinusoidal on the second
%       row, and P vs.E on the third row
%   varmat: a structure with the variance-covariance matrix from the three
%       models. The variance-covariance matrix is obtained using the
%       parameters obtained from bootstrapping the residuals
%   fitline: a table with the model line obtained using the three models
%       and the vector of the dates with 240 points per day
% EXAMPLE:
%   x = 736497 + [0.1 0.25 0.4 0.7 1.1 1.25 1.35 1.65 1.8 2.05 2.15 2.45 2.6 2.9];
%   y = 200 + [0.121 -0.240 -0.223 0.259 -0.017 -0.051 -0.107 0.307 0.324 0.063 -0.026 -0.226 0.426 0.373];
%   [rates,residuals,varmat,fitline] = dielFit(x,y,22.75,-158,-10);
%
% ACKNOWLEDGEMENT:
% This routine builds on code provided by David P. Nicholson (WHOI).
%
% REFERENCES:
% Draper, N. R., and H. Smith, (2014), Applied regression analysis, Third
% ed. John Wiley & Sons. doi: 10.1002/9781118625590
%
% Durbin, J., and G. S. Watson, (1950), Testing for serial correlation in
% least squares regression. I, Biometrika, 37, 409?428,
% doi: 10.1093/biomet/37.3-4.409.
%
% Nicholson, D. P., S. T. Wilson, S. C. Doney, and D. M. Karl (2015),
% Quantifying subtropical North Pacific gyre mixed layer primary
% productivity from Seaglider observations of diel oxygen cycles. Geophys.
% Res. Lett., 42, 4032?4039. doi: 10.1002/2015GL063065.
%
% VERSION HISTORY
%   r0.11: First release (non-linear least squares)
%   r0.12: Adds computation of 95% confidence intervals using nlparci.m
%   r0.13: The percent level for confidence intervals can now be assigned
%   r0.14: Change method to linear least squares with nonnegativity constraints
%          Bootstrapping residuals to compute confidence intervals
%          Adds p-value from Durbin-Watson test (residual autocorrelation)
%   r0.15: Change method to ordinary linear least squares
%          Input dates and conc now accept NaNs that are removed from the
%           vector before the fit
%   r0.16: Output 'rates' now contains third fit parameter, C0
%          New output variable 'fitline' with lines from fit results
%   r0.17: Computation of variance-covariance matrix (output 'varmat')
%          'alpha' default is now 0.317
%   r0.18: Computation of variance-covariance matrix is now done by
%          bootstrapping the residuals
%
% Benedetto Barone - September 2018 - Revision 0.18

p = inputParser;

defaultAlpha = 0.317; % 68.3% confidence intervals
defaultEmaxEk = 1; % no photoinhibition
defaultEmaxEb = 0; % Ek == maximum daily irradiance

defaultMethod = 'sin';
validMethod = {'linear','sin','par'};
checkMethod = @(x) any(validatestring(x,validMethod));


addRequired(p,'mtime',@isnumeric);
addRequired(p,'conc',@isnumeric);
addRequired(p,'lat',@isnumeric);
addRequired(p,'lon',@isnumeric);
addParameter(p,'alpha',defaultAlpha,@isnumeric);
addParameter(p,'EmaxEk',defaultEmaxEk,@isnumeric);
addParameter(p,'EmaxEb',defaultEmaxEb,@isnumeric);
addParameter(p,'method',defaultMethod,checkMethod);


parse(p,mtime,conc,lat,lon,varargin{:});

alpha = p.Results.alpha;
EmaxEk = p.Results.EmaxEk;
EmaxEb = p.Results.EmaxEb;
method = p.Results.method;

% get rid of any nans
gd = ~isnan(mtime + conc);
xfit = mtime(gd) - fix(min(mtime(gd))); yfit = conc(gd);

nt = 240; % number of time points per day
% Extract solar elevation cycle
timezone = 0;
[~,t,~,z] = suncycle(lat,lon,mean(mtime,'omitnan') - timezone/24,nt);
% Transform t 
td = t - floor(min(t));
ndays = ceil(max(td))+1;
%[t,ind_row] = sort(t);
% Normalized light intensity from solar elevation (1 is maximum light, Emax)
%z = z(ind_row);
z(z < 0) = 0;
Erel = sind(z);

% edges of tbins
tt = linspace(0,ndays,nt*ndays+1)';



% 1. Linear production model
if strcmpi(method,'linear')
    P = zeros*td;
    P(Erel>0) = 1;
    
    % 2. Sinusoidal producion model (linear with light)
elseif strcmpi(method,'sin')
    P = Erel;
elseif strcmpi(method,'par')
    % 3. PvsE producion model
    P = (1-exp(-EmaxEk.*Erel)).*exp(-EmaxEb.*Erel);
else
    error('oops')
end
Pt = interp1(td,P,tt,'linear','extrap');
Ptn = Pt./trapz(tt,Pt);
A = [ones(length(xfit),1) interp1(tt,cumtrapz(tt,Ptn),xfit)' -xfit'];

% Fit moddels
par = A\yfit';
% Residuals
resid = yfit'-A*par;
% Y predicted from the model fit
y_hat = A*par;

% Bootstrap residual to compute parameter confidence intervals (200 iterations)
boot = bootstrp(200,@(bootr) A\(y_hat+bootr),resid);

% Confidence intervals
ci = [prctile(boot,(alpha/2)*100)' prctile(boot,(1-alpha/2)*100)'];

% Model statistics (R^2 and p value)
[r_temp,p_temp] = corrcoef(yfit,y_hat);
rsq = r_temp(2)^2;
pval = p_temp(2);

% Variance-covariance matrix from bootstrap
varmat = cov(boot);
% Durbin-Watson test for residual autocorrelation
pdw = dwtest(resid,A);
% Output variables
rates = table;
rates.GPP = par(2);
rates.CR = par(3);
rates.C0 = par(1);
rates.R2 = rsq;
rates.p = pval;
rates.pdw = pdw;
residuals = resid;


end