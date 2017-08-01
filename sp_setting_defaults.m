function spParams = sp_setting_defaults(sp,spParamsUser,srate)
% Establish basic parameter settings, then update for species
% specific defaults, and user preferences.
% Pulled into subroutine kf 10/4/2016

spParams = [];

specChar = 'Unk';  %Simone abbreviation for species
speName = 'Unknown';  % Species code used in file names 
tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
dtHi = .5; % max yaxis value for IPI display in sec
fLow = 5; % boundary for spectrum plot
threshRL = 0; % minimum RL threshold in dB peak-to-peak
threshRMS = 0; % default for < coommand, RMS threshold cutoff
threshPP = 0; % default for : coommand, PP threshold cutoff
threshHiFreq = 0; % default for ^ command, high freq cutoff for clicks
ltsaContrast = 250; % ltsa contrast
ltsaBright = 100; % ltsa brightness
ltsaLims = [0,100]; % max and min of LTSA plot
ltsaMax = 6; % ltsa maximum duration per session
rlLow = 110; % PP plot window low limit
rlHi = 170; % PP plot window high limit
dfManual = []; % LTSA step size in 10 [Hz] bins
p1Low = threshRL - 5;
p1Hi = 170; % ??
minSpectralFreq = 0; 
minBout = [];
slope = 1;


if (strcmp(sp,'Ko') || strcmp(sp,'k'))
    specChar = 'K';
    speName = 'Kogia';
    tfSelect = 80000; 
    dtHi = 0.5;  
    fLow = 70;   
    threshRL = 116;
elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
    specChar = 'Z';
    speName = 'Cuviers';  tfSelect = 40200;
    dtHi = 1.0; 
    fLow = 25;  
    threshRL = 121; 
    ltsaContrast = 200; ltsaBright = 30;
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    specChar = 'M';
    speName = 'Gervais';
    tfSelect = 40200;
    dtHi = 1.0; 
    fLow = 25; 
    threshRL = 121;
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    specChar = 'G'; %Simone abbreviations for BW species
    speName = 'BWG';  tfSelect = 40200; % freq used for transfer function
    dtHi = 1.0; % scale for IPI display in sec
    fLow = 25;   % 25 kHz boundary for spec plot
    threshRL = 121; % dB threshold
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    specChar = 'D'; 
    speName = 'BW31';  tfSelect = 40200; 
    dtHi = 1.0;
    fLow = 25; 
    threshRL = 121;
elseif strcmpi(sp,'De')
    speName = 'Delphin';  
    tfSelect = 0; % already in dB no correction
    dtHi = 0.6; 
    fLow = 10;   
    threshRL = 118;
    rlLow = threshRL - 6.9; rlHi = 190;
elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
    speName = 'Porpoise';  tfSelect = 0; 
    dtHi = 0.5; 
    fLow = 25;  
    threshRL = 100;
    rlLow = threshRL - 5; rlHi = 190;
    ltsaContrast = 250; ltsaBright = 100; 
elseif strcmpi(sp,'MFA')
    speName = 'MFA';  tfSelect = 4000; 
    ltsaMax = .5; 
    threshRL = 80;
    dtHi = 2;
    fLow = 2;  
    rlLow = threshRL - 5; rlHi = 180;
    dfManual = 10;   
elseif strcmpi(sp,'whs')
    speName = 'whs';  tfSelect = 0; 
    ltsaMax = .5;
    threshRL = 0;
    dtHi = 2; 
    fLow = 5;   
    rlLow = 0; rlHi = 20;
    dfManual = 10;   
    ltsaContrast = 310; ltsaBright = 100; 
    ltsaLims = [5,30];
elseif strcmpi(sp,'Dl')
    speName = 'Beluga'; tfSelect = 45000;
    dtHi = 0.5;
    fLow = 20;  
    threshRL = 110;
    rlLow = threshRL - 5; rlHi = 170;
    ltsaContrast = 200; ltsaBright = 70;
elseif (strcmp(sp,'PM') || strcmp(sp,'pm') ...
        || strcmp(sp,'Pm'))
    speName = 'Pm'; tfSelect = 15000;
    dtHi = 2; 
    fLow = 5;
    threshRL = 120; threshHiFreq = 20;
    threshRMS = 95; threshPP = 145;
    ltsaContrast = 180; ltsaBright = 73;
    dfManual = 100;
    minBout = 75;
    slope = 1.2;
else
    warning('Unknown Species Type!!!')
end

%% apply default if user has not specified a value
if isfield(spParamsUser,'specChar')
    spParams.specchar = spParamsUser.specChar;
else
    spParams.specchar = specChar;
end

if isfield(spParamsUser,'speName')
    spParams.speName = spParamsUser.speName;
else
    spParams.speName = speName;
end

if isfield(spParamsUser,'tfSelect')
    spParams.tfSelect = spParamsUser.tfSelect;
else
    spParams.tfSelect = tfSelect;
end

if isfield(spParamsUser,'dtHi')
    spParams.dtHi = spParamsUser.dtHi;
else
    spParams.dtHi = dtHi;
end

if isfield(spParamsUser,'fLow')
    spParams.fLow = spParamsUser.fLow;
else
    spParams.fLow = fLow;
end

if isfield(spParamsUser,'threshRL')
    spParams.threshRL = spParamsUser.threshRL; 
else
    spParams.threshRL = threshRL;
end

if isfield(spParamsUser,'threshRMS')
    spParams.threshRMS = spParamsUser.threshRMS; 
else
    spParams.threshRMS = threshRMS;
end

if isfield(spParamsUser,'threshPP')
    spParams.threshPP = spParamsUser.threshPP; 
else
    spParams.threshPP = threshPP;
end

if isfield(spParamsUser,'threshHiFreq')
    spParams.threshHiFreq = spParamsUser.threshHiFreq; 
else
    spParams.threshHiFreq = threshHiFreq;
end

if isfield(spParamsUser,'ltsaContrast')
    spParams.ltsaContrast = spParamsUser.ltsaContrast;    
else
    spParams.ltsaContrast = ltsaContrast;
end

if isfield(spParamsUser,'ltsaBright')
    spParams.ltsaBright = spParamsUser.ltsaBright;
else
    spParams.ltsaBright = ltsaBright;
end

if isfield(spParamsUser,'ltsaLims')
    spParams.ltsaLims = spParamsUser.ltsaLims;
else
    spParams.ltsaLims = ltsaLims;
end

if isfield(spParamsUser,'ltsaMax')
    spParams.ltsaMax = spParamsUser.ltsaMax;    
else
    spParams.ltsaMax = ltsaMax;
end

if isfield(spParamsUser,'rlLow')
    spParams.rlLow = spParamsUser.rlLow;
else
    spParams.rlLow = rlLow;
end

if isfield(spParamsUser,'rlHi')
    spParams.rlHi = spParamsUser.rlHi;  
else
    spParams.rlHi = rlHi;
end

if isfield(spParamsUser,'dfManual')
    spParams.dfManual = spParamsUser.dfManual;    
else
    spParams.dfManual = dfManual;
end

if isfield(spParamsUser,'dfManual')
    spParams.dfManual = spParamsUser.dfManual;
else
    spParams.dfManual = dfManual;
end

if isfield(spParamsUser,'p1Low')
    spParams.p1Low = spParamsUser.p1Low;
else
    spParams.p1Low = p1Low;
end

if isfield(spParamsUser,'p1Hi')
    spParams.p1Hi = spParamsUser.p1Hi;
else
    spParams.p1Hi = p1Hi;
end

if isfield(spParamsUser,'minBout')
    spParams.minBout = spParamsUser.minBout;
else
    spParams.minBout = minBout;
end

if isfield(spParamsUser,'slope')
    spParams.slope = spParamsUser.slope;
else
    spParams.slope = slope;
end
if isfield(spParamsUser,'minSpectralFreq')
    spParams.minSpectralFreq = spParamsUser.minSpectralFreq;
else
    spParams.minSpectralFreq = minSpectralFreq;
end

if isfield(spParamsUser,'maxSpectralFreq')
    spParams.maxSpectralFreq = spParamsUser.maxSpectralFreq;
else
    spParams.maxSpectralFreq = srate/2; % set it later based on sample rate.
end

if isfield(spParamsUser,'minSpectralFreq')
    spParams.minSpectralFreq = spParamsUser.minSpectralFreq;
else
    spParams.minSpectralFreq = minSpectralFreq;
end

if isfield(spParamsUser,'maxSpectralFreq')
    spParams.maxSpectralFreq = spParamsUser.maxSpectralFreq;
else
    spParams.maxSpectralFreq = srate/2; % set it later based on sample rate.
end
