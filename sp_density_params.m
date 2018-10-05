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
    if strcmp(site,'GC')
        fpRate = 0.0274;
        clickRate = 1.4454;
        pDet = 0.0356;
        maxRadius_km = 12;
    elseif strcmp(site,'MC980')
        fpRate = 0.0649;
        clickRate = 1.4982;
        pDet = 0.06069;
        maxRadius_km = 12;
    elseif strcmp(site,'MC800')
        fpRate = 0.0649;
        clickRate = 1.4982;
        pDet = 0.06384;
        maxRadius_km = 12;
    elseif strcmp(site,'DT')
        fpRate = 0.0815;
        clickRate = 1.3515;
        pDet = 0.05276;
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
spParams.clickRate = clickRate;
spParams.pDet = pDet;
spParams.maxRadius_km = maxRadius_km;
