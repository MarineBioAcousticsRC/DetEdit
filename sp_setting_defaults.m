function spParams = sp_setting_defaults(sp,srate,analysis)
% Establish basic parameter settings, then update for species
% specific defaults, and user preferences.
% Pulled into subroutine kf 10/4/2016

% Set default parameters
spParams = [];
specChar = 'Unk';  %Simone abbreviation for species
speName = 'Unknown';  % Species code used in file names 
tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
dtHi = .5; % max yaxis value for IPI display in sec
fLow = 0; % Minimum frequency of interest
fHi = srate/2; % Maximum frequency of interest
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
minBout = [];% minimum bout duration
gth = .5;    % gap time in hrs between sessions
minDur = []; % minimum window duration (if specified in minutes)
slope = 1; % slope for shifting data in plots 51 and 53
% all parameters for modDet
iciRange = []; % min/max ici in ms for modDet plots
dbRange = [];  % min/max db for modDet plots of pp and rms
frRange = [fLow, fHi];   % min/max frequency for modDet plots of peak and center freq
frdbwRange = [fLow, fHi]; % min/max frequency for modDet plots of 3/10 db bw
durRange = []; % min/max duration in us for modDet plots
durstep = 1; % step range for number bins in histogram 
N = srate; % FFT size for parameter clicks

% Set parameters according to sp
if (strcmp(sp,'Ko') || strcmp(sp,'k'))
    speName = 'Kogia'; specChar = 'K'; 
    tfSelect = 80000; 
    dtHi = 0.5;  
    fLow = 70;   
    threshRL = 116;
    iciRange = [40, 130];
    dbRange = [100, 150];
    frRange = [80, 160];
    frdbwRange = [0, 80];
    durRange = [30, 111];
    durstep = 3;
elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
    speName = 'Cuviers'; specChar = 'Z'; 
    tfSelect = 40200;
    dtHi = 1.0; 
    fLow = 25;  
    threshRL = 121; 
    ltsaContrast = 200; ltsaBright = 30;
    iciRange = [40, 750];
    dbRange = [90, 170];
    frdbwRange = [0, 80];
    durRange = [30, 300];
    durstep = 2;
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    speName = 'Gervais'; specChar = 'M'; 
    tfSelect = 40200;
    dtHi = 1.0; 
    fLow = 25; 
    threshRL = 121;
    iciRange = [40, 750];
    dbRange = [90, 170];
    frdbwRange = [0, 80];
    durRange = [30, 300];
    durstep = 2;
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    speName = 'BWG'; specChar = 'G'; 
    tfSelect = 40200;
    dtHi = 1.0; 
    fLow = 25;   
    threshRL = 121; 
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    speName = 'BW31'; specChar = 'D'; 
    tfSelect = 40200; 
    dtHi = 1.0;
    fLow = 25; 
    threshRL = 121;
elseif strcmpi(sp,'De')
    speName = 'Delphin';  
    dtHi = 0.6; 
    fLow = 10;   
    threshRL = 118;
    rlLow = threshRL - 6.9; rlHi = 190;
elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
    speName = 'Porpoise'; 
    dtHi = 0.5; 
    fLow = 25;  
    threshRL = 100;
    rlLow = threshRL - 5; rlHi = 190;
    ltsaContrast = 250; ltsaBright = 100; 
elseif strcmpi(sp,'MFA')
    speName = 'MFA';  
    tfSelect = 4000; 
    ltsaMax = .5; 
    threshRL = 80;
    dtHi = 2;
    fLow = 2;  
    rlLow = threshRL - 5; rlHi = 180;
    dfManual = 10;   
elseif strcmpi(sp,'whs')
    speName = 'whs';
    ltsaMax = .5;
    dtHi = 2; 
    fLow = 5;   
    rlLow = 0; rlHi = 20;
    dfManual = 10;   
    ltsaContrast = 310; ltsaBright = 100; 
    ltsaLims = [5,30];
elseif strcmpi(sp,'Dl')
    speName = 'Beluga'; 
    tfSelect = 45000;
    dtHi = 0.5;
    fLow = 20;  
    threshRL = 110;
    rlLow = threshRL - 5; rlHi = 170;
    ltsaContrast = 200; ltsaBright = 70;
    iciRange = [20, 500];
    dbRange = [90, 170];
    durRange = [10, 300];
    durstep = 2;
elseif (strcmp(sp,'PM') || strcmp(sp,'pm') || strcmp(sp,'Pm'))
    speName = 'Pm'; 
    dtHi = 2; 
    fLow = 5;
    threshRL = 120; threshHiFreq = 20;
    threshRMS = 95; threshPP = 145;
    ltsaContrast = 180; ltsaBright = 73;
    dfManual = 100;
    minBout = 75; minDur = 60; 
    slope = 1.2;
    iciRange = [50, 4000];
    frRange = [fLow fHi];
    N = 512;
else
    warning('Unknown Species Type!!!')
end

%% create struct to return parameters
switch analysis
    case {'detEdit','mkLTSA'} 
        spParams.specchar = specChar;
        spParams.speName = speName;
        spParams.tfSelect = tfSelect;
        spParams.dtHi = dtHi;
        spParams.fLow = fLow;
        spParams.fHi = fHi;
        spParams.threshRL = threshRL;
        spParams.threshRMS = threshRMS;
        spParams.threshPP = threshPP;
        spParams.threshHiFreq = threshHiFreq;
        spParams.ltsaContrast = ltsaContrast;
        spParams.ltsaBright = ltsaBright;
        spParams.ltsaLims = ltsaLims;
        spParams.ltsaMax = ltsaMax;
        spParams.rlLow = rlLow;
        spParams.rlHi = rlHi;
        spParams.dfManual = dfManual;
        spParams.dfManual = dfManual;
        spParams.p1Low = p1Low;
        spParams.p1Hi = p1Hi;
        spParams.minBout = minBout;
        spParams.minDur = minDur;
        spParams.slope = slope;
    case {'modDet'} 
        spParams.tfSelect = tfSelect;
        spParams.dbRange = dbRange;
        spParams.iciRange = iciRange;
        spParams.frRange = frRange;
        spParams.frdbwRange = frdbwRange;
        spParams.durRange = durRange;
        spParams.durstep = durstep;
        spParams.N = N;
        spParams.gth = gth;
        spParams.minBout = minBout;
    otherwise
        sprintf(['No analysis specified. Please add one of these options:\n',...
        'detEdit, mkLTSA or modDet'])
end



