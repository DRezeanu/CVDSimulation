function sensitivity = neitzSpec(peak, pigmentOD, varargin)

%NEITZSPEC - Photopigment template function developed in the Neitz Lab at
%University of Washington
%   neitzSpec(X, Y) computes the log spectral sensitivity for a photopigment
%   with a given peak wavelength X and an optical density Y, evaluated
%   between 400 and 700nm, with a step of 1nm. 
%
%   neitzSpec(X,Y,Z) computes the log spectral sensitivity for a photopigment
%   with a given peak wavelength X and an optical density Y, evaluated at
%   the wavelengths defined by the vector Z. NOTE: Z must be between
%   390 and 800, with a minimum step of 0.1nm.
%
%   neitzSpec(..., SCALE) computes the spectral sensitivity using the
%   scaling defined by SCALE:
%       'log' - (default) returns log10 spectral sensitivity curve
%       'linear' - returns linear spectral sensitivity curve
%
%   neitzSpec(..., UNITS) defines whether the plot will be in units of
%   energy or quanta.
%       'quanta' - (default) returns spectral sensitivity in quanta
%       'energy' - returns spectral sensitivity in energy
%
%   neitzSpec(..., AGE) scalar value defines the age used to compute
%   occular media optical density spectra. Only used if lensOD is set 
%   to 'vanNorren'. Default value is 35.
%
%   neitzSpec(..., LENSOD) computes spectral sensitivity using the lens
%   optical density filtering defined by LENSOD:
%       'Stockman' - (default) returns spectral sensitivity filtered using
%       the lens optical density spectrum available from CVRL.org
%       'vanNorren' - returns spectral sensitivity filtered using age-based 
%       lens and ocular media density spectrum computed using van Norren 
%       and Van de Kraats equation.
%
%   neitzSpec(..., MACULAOD) scalar value between 0 and 1 that defines the
%   proportion of macular pigment filtering to be used when applying the
%   Stockman macular pigment optical density spectrum from CVRL.org.
%   Default value is 1, representing full macular pigment filtering in the
%   very center of the fovea. Lower values correspond to decreased macular
%   pigment filtering with increased eccentricity.
%
%   neitzSpec(..., PRERETINALFILTERS) logical value, true by default. If
%   false, the function will not account for the filtering of the lens and
%   macular pigment.
%
%   neitzSpec(..., CHECKSPECTRUM) logical value, false by default. If true,
%   the function will plot the output spectrum.

        defaultWavelength = 400:1:700;
        defaultAge = 0;
        defaultMaculaOD = 1;
        defaultScale = 'linear';
        expectedScales = {'log', 'linear'};
        defaultLensOD = 'Stockman';
        expectedLenses = {'Stockman', 'vanNorren'};
        defaultUnits = 'quanta';
        expectedType = {'quanta', 'energy'};
        defaultPreRetinalFilters = true;
        defaultCheckSpectrum = false;
        
        ip = inputParser;

        validPeak = @(x) isscalar(x) && (x >= 390) && (x <= 800);
        validPigmentOD = @(x) isscalar(x) && (x >= 0) && (x <= 1);
        validWavelengths = @(x) isvector(x) && (min(x) >= 390) && (max(x) <=800);
        validAge = @(x) isscalar(x) && (x>0);
        validMaculaOD = @(x) (x>0) && (x<=1);
        validLensOD = @(x) any(validatestring(x, expectedLenses));
        validScale = @(x) any(validatestring(x,expectedScales));
        validUnits = @(x) any(validatestring(x, expectedType));
        
        addRequired(ip, 'peak', validPeak);
        addRequired(ip, 'pigmentOD', validPigmentOD);
        addOptional(ip, 'wavelength', defaultWavelength, validWavelengths);
        addParameter(ip, 'age', defaultAge, validAge);
        addParameter(ip, 'maculaOD', defaultMaculaOD, validMaculaOD);
        addParameter(ip, 'lensOD', defaultLensOD, validLensOD);
        addParameter(ip, 'scale', defaultScale, validScale);
        addParameter(ip, 'units', defaultUnits, validUnits);
        addParameter(ip, 'checkSpectrum', defaultCheckSpectrum, @(x) islogical(x));
        addParameter(ip, 'preRetinalFilters', defaultPreRetinalFilters, @(x) islogical(x));
        
        parse(ip, peak, pigmentOD, varargin{:})
        

        % Pull Validated Values from Input Parser
        
        wavelength = ip.Results.wavelength;
        age = ip.Results.age;
        maculaOD = ip.Results.maculaOD;
        lensOD = ip.Results.lensOD;
        scale = ip.Results.scale;
        checkSpectrum = ip.Results.checkSpectrum;
        units = ip.Results.units;
        preRetinalFilters = ip.Results.preRetinalFilters;
        

        % Extinction Contstants
        
        peakSensitivity=(log10(1/peak)-log10(1/558.5));
        constants = readmatrix('Photopigment Template.xlsx', 'Range', 'A1:Y1');
        constants = [constants, peakSensitivity];
        

        % Load or Derive Lens and Macular Pigment Absorption
        
        stockmanLens_full = readmatrix('stockmanLens.csv', 'Range', 'A1:B4401');
        stockmanLens = zeros(1,length(wavelength));
        
        for i = 1:length(wavelength)
            stockmanLens(i) = stockmanLens_full(find(stockmanLens_full(:,1) == wavelength(i)),2);
        end
        
        if contains(lensOD, 'Stockman', 'IgnoreCase', true)
            lens = stockmanLens;
        else    
            lens = lensAge(age, wavelength);
        end
        
        stockmanMacula_full = readmatrix('stockmanMacula.csv', 'Range', 'A1:B4401');
        stockmanMacula = zeros(1, length(wavelength));
        
        for i = 1:length(wavelength)
            stockmanMacula(i) = stockmanMacula_full(find(stockmanMacula_full(:,1) == wavelength(i)), 2);
        end
        
        macula = stockmanMacula.*maculaOD;
        

        %  Calculate Extinction 
        
        extinction=zeros(1,length(wavelength));
        
        for x=1:length(wavelength)
            extinction(x)=log10(-constants(5)+constants(5)*tanh(-(((10^(log10(1/wavelength(x))...
                -constants(26))))-constants(6))/constants(7)))+constants(4)...
                +constants(1)*tanh(-(((10^(log10(1/wavelength(x))-constants(26))))...
                -constants(2))/constants(3))-(constants(10)/constants(9)...
                *((1/(sqrt(2*pi)))*exp(1)^(-0.5*(((10^(log10(1/wavelength(x))...
                -constants(26)))-constants(8))/constants(9))^2)))...
                -(constants(13)/constants(12)*((1/(sqrt(2*pi)))...
                *exp(1)^(-0.5*(((10^(log10(1/wavelength(x))-constants(26)))...
                -constants(11))/constants(12))^2)))-(constants(16)/constants(15)...
                *((1/(sqrt(2*pi)))*exp(1)^(-0.5*(((10^(log10(1/wavelength(x))...
                -constants(26)))-constants(14))/constants(15))^2)))+(constants(19)...
                /constants(18)*((1/(sqrt(2*pi)))*exp(1)^(-0.5*(((10^(log10(1/wavelength(x))...
                -constants(26)))-constants(17))/constants(18))^2)))...
                +((constants(22)/constants(21)*((1/(sqrt(2*pi)))*exp(1)^(-0.5...
                *(((10^(log10(1/wavelength(x))-constants(26)))-constants(20))/constants(21))...
                ^2)))/10)+((constants(25)/constants(24)*((1/(sqrt(2*pi)))...
                *exp(1)^(-0.5*(((10^(log10(1/wavelength(x))-constants(26)))...
                -constants(23))/constants(24))^2)))/100);
        end
        

        % Account for Optical Density of Photopigment
 
        extinction_od = log10((1-10.^-((10.^extinction).*pigmentOD))./(1-10.^-...
            pigmentOD));
        
        % Convert to energy if plotting energy

        if contains(units, 'energy', 'IgnoreCase', true)
            conversionFactor = wavelength./wavelength(end);
            extinction_od = extinction_od+log10(conversionFactor);
        end

        % Account for preretinal filters if "preRetinalFilters" is set to true

        if preRetinalFilters == true
            extinction_od = extinction_od-lens-macula;
        end

        % Calculate Spectral Sensitivity
        
        if contains(scale, 'log', 'IgnoreCase', true)
            sensitivity = extinction_od;
        else
            sensitivity=10.^(extinction_od);
        end


        % Plot output if "checkSpectrum" is set to true

        if checkSpectrum == true
            specPlot = figure;
            specPlot.Position = [0,0, 1920/2, 1080/2];

            plot(wavelength, sensitivity, '-k', 'LineWidth', 2);
            grid on
            xlabel('Wavelength (nm)')
            if contains(scale, 'log', 'IgnoreCase', true)
                ylabel('Log Sensitivity')
            else
                ylabel('Sensitivity')
            end
        end

end