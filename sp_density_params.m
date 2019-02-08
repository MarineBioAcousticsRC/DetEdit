function spParams = sp_density_params(varargin)
% Establish density parameter, then update for species
% specific defaults, and user preferences.
% Pulled into subroutine asb 10/5/2018

% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'sp'
            sp = varargin{n+1}; n=n+2;
        case 'site'
            site = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

% Set default parameters
fpRate = 0; % False Positive Rate
clickRate = []; % Click Rate (clicks/second)
pDet = []; % Probability of detection
maxRadius_km = []; % Maxium range (km)

%% Set parameters according to sp
% sperm whale parameters for density
if (strcmp(sp,'PM') || strcmp(sp,'pm') || strcmp(sp,'Pm'))
    %     if strcmp(site,'GC')
    %         fpRate = 0.0274;
    %         fpRateCV = 0.2255;
    %         clickRate = 1.117;
    %         clickRateCV = 0.4029;
    %         pDet = 0.0554;
    %         pDetCV = 0.4880;
    %         maxRadius_km = 12;
    if strcmp(site,'GC')
        fpRate = 0.0274;
        fpRateCV = 0.2255;
        clickRate = 1.4454;
        clickRateCV = 0.3478;
        pDet = 0.0558;
        pDetCV = 0.4836;
        maxRadius_km = 12;
    elseif strcmp(site,'MC980')
        fpRate = 0.0669;
        fpRateCV = 0.1409;
        clickRate = 1.5129;
        clickRateCV = 0.3478;
        pDet = 0.0648;
        pDetCV = 0.4786;
        maxRadius_km = 12;
    elseif strcmp(site,'MC800')
        fpRate = 0.0343;
        fpRateCV = 0.1968;
        clickRate = 1.4493;
        clickRateCV = 0.3478;
        pDet = 0.0632;
        pDetCV = 0.4629;
        maxRadius_km = 12;
    elseif strcmp(site,'DTnoSeas')
        fpRate = 0.0815;
        fpRateCV = 0.5198;
        clickRate = 1.3515;
        clickRateCV = 0.3478;
        pDet = 0.0533;
        pDetCV = 0.5007;
        maxRadius_km = 12;
    elseif strcmp(site,'DTinSeas')
        fpRate = 0.0815;
        fpRateCV = 0.5198;
        clickRate = 1.3671;
        clickRateCV = 0.3478;
        pDet = 0.0531;
        pDetCV = 0.4943;
        maxRadius_km = 12;
    elseif strcmp(site,'DToutSeas')
        fpRate = 0.0815;
        fpRateCV = 0.5198;
        clickRate = 1.2824;
        clickRateCV = 0.3478;
        pDet = 0.0518;
        pDetCV = 0.4857;
        maxRadius_km = 12;
    else
        warning('Unknown Site!!!')
    end
    
    %%%% "Species" parameters for density goes here
    
else
    warning('Unknown Species Type!!!')
end

%% create struct to return parameters
spParams.fpRate = fpRate;
spParams.fpRateCV = fpRateCV;
spParams.clickRate = clickRate;
spParams.clickRateCV = clickRateCV;
spParams.pDet = pDet;
spParams.pDetCV = pDetCV;
spParams.maxRadius_km = maxRadius_km;
