function lensOD = lensAge(age, varargin)
%LENSAGE - Age-based lens and ocular media optical density spectrum
%   lensAge(X) computes the lens and ocular media optical density spectrum 
%   for a person of age X, evaluated between 400 and 700nm, with a step of
%   1nm, using the van Norren and Van de Kraats equation.
%
%   lensAge(X,Y) computes the lens and ocular media optical density spectrum
%   for a person of age X, evaluated between the values defined by the
%   vector Y. NOTE: Y must be between 390 and 800nm, with a minimum step of
%   0.1nm.
%
%   lensAge(..., CHECKSPECTRUM) logical value, default is false. If true,
%   the function will plot the output spectrum


defaultWavelength = 400:1:700;
defaultCheckSpectrum = false;

validAge = @(x) isscalar(x) && (x>0) && (x < 100);
validWavelength = @(x) isvector(x) && (min(x) >= 390) && (max(x) <= 800);

ip = inputParser;

addRequired(ip, 'age', validAge);
addOptional(ip, 'wavelength', defaultWavelength, validWavelength);
addParameter(ip, 'checkSpectrum', defaultCheckSpectrum, @(x) islogical(x))

parse(ip, age, varargin{:})

WL = ip.Results.wavelength;
checkSpectrum = ip.Results.checkSpectrum;

densityConstants = [0.446, 0.000031, 14.19, 10.68, 0.057, 273, 0.998, 0.000063, 2.13, 0.029, 370, 0.059, 0.000186, 11.95, 0.021, 325, 0.016, 0.000132, 1.43, 0.008, 325, 0.111];

lensOD = zeros(1,length(WL));

for x = 1:length(WL)
    lensOD(x) = (densityConstants(1)+densityConstants(2)*age^2)*(400/WL(x))^4+densityConstants(3)*densityConstants(4)*exp(-((densityConstants(5)*(WL(x)-densityConstants(6)))^2))+...
        (densityConstants(7)-densityConstants(8)*age^2)*densityConstants(9)*exp(-((densityConstants(10)*(WL(x)-densityConstants(11)))^2))+(densityConstants(12)+densityConstants(13)*age^2)...
        *densityConstants(14)*exp(-((densityConstants(15)*(WL(x)-densityConstants(16)))^2))+(densityConstants(17)+densityConstants(18)*age^2)*densityConstants(19)...
        *exp(-((densityConstants(20)*(WL(x)-densityConstants(21)))^2))+densityConstants(22);
end

if checkSpectrum == true
    lensPlot = figure;
    lensPlot.Position = [0,0, 720, 720];

    plot(WL, lensOD, '--k', 'LineWidth', 2)
    grid on
    xlabel('Wavelength (nm)')
    ylabel('Log Optical Density')
end