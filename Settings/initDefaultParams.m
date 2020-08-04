%initDefaultParams

% Initialize default string parameters for the interface, mkLTSAsessions,
% modDet and summaryParams

sampleRate = 200;
speName = '';
ltsaDir = '';
tfName = '';
effortTimes = '';
referenceTime = '';

% General bout parameters

threshRL = 0; % minimum RL threshold in dB peak-to-peak
tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
minBout = 0;% minimum bout duration in seconds
gth = .5;    % gap time in hrs between sessions
binDur = 5; % bin duration in minutes
dfManual = []; % LTSA step size in 10 [Hz] bins
specploton = 1; %(1 = yes | 0 = no ) spectra parameters plot
minNdet = 1; % minimum number of detections per session. Sessions with fewer than this will be skipped
maxDetLoad = 4e5; % [] - read all or 4e5 - the number of detections above 
% which you want to read from disk instead of loading all spectra and 
% timeseries into memory this is for large files (e.g. dolphin click detections)
% if maxDetLoad exist, plotaxes can be defined to keep the format of the
% panels        
c4fd = 1; % Detections step size to estimate false detection rate
rawFileDur = []; % Raw file length, default 75s


% Parameters for the interface

% Panel LTSA and time series
rlLow = 110; % PP plot window low limit
rlHi = 170; % PP plot window high limit
ltsaContrast = 250; % ltsa contrast
ltsaBright = 100; % ltsa brightness
ltsaLims = [0,sampleRate/2]; % max and min of LTSA plot
ltsaMax = 6; % ltsa maximum duration per session
dtHi = .5; % max yaxis value for ICI display in sec
minDur = []; % minimum window duration (if specified in minutes)

% Panel Frequency spectra
fLow = 0; % Minimum frequency of interest
fHi = sampleRate/2; % Maximum frequency of interest

% Panel RL rms vs. RL pp | Peak freq.
slope = 1; % slope for shifting data vertically
rmsLow = 80; % transformed received level (dBrms) plot window low limit
rmsHi = 140; % transformed received level (dBrms) plot window high limit
threshRMS = 0; % default for < command, RMS threshold cutoff
threshPP = 0; % default for : command, PP threshold cutoff
threshHiFreq = 0; % default for ^ command, high freq cutoff for clicks
autoFalse = false; % Apply automatic false thresholds to entire file. 
 
% Colors to use for classification - ID signal types
colorTab = [204, 204, 255; ... % Blainville's, lilac
            255, 0, 0;... % Boats, red
            255, 128, 0; ... % CT11, orange
            102, 255, 178; ... % CT2+CT9, seafoam
            0, 153, 0; ... % CT3+CT7, crayola green
            0, 128, 255; ... % CT4/6+CT10, medium blue
            51, 255, 51; ... % CT5, bright green
            102, 178, 255; ... % CT8, periwinkle
            0, 255, 255; ... % Cuvier's, cyan
            255, 204, 153; ... % Gervais, tan
            245, 194, 66;... % GoM Gervais, mustard yellow
            255, 0, 0;... % HFA, red
            153, 51, 255; ... % Kogia, purple
            255, 0, 0;... % MFA, red
            255, 0, 0;... % MultiFreq_Sonar, red
            204, 255, 153; ... % Rissos, lime green
            255, 0, 0;... % Snapping Shrimp, red
            255, 0, 255; ... % Sowerby's, magenta
            102, 51, 0; ... % Sperm whale, brown
            249, 177, 211]./255; % True's, light pink
colorTab = round(colorTab.*100)/100;      

mySpID = {'Blainvilles','Boats','CT11','CT2+CT9','CT3+CT7','CT4/6+CT10',...
    'CT5','CT8','Cuviers','Gervais','GoM_Gervais','HFA','Kogia','MFA',...
    'MultiFreq_Sonar','Rissos','Snapping_Shrimp','Sowerbys','Sperm_Whale',...
    'Trues'};
minLabelConfidence = 0; % minimum classifciation confidence threshold for 
% displaying signal labels


%% Parameters for modDet (if specified to compute)

excludeID = 0; % exclude ID times from MTT files
calcParams = 0; % yes - 1 | no - 0. Calculate Parameters peak-to-peak, 
% inter-detection-interval and peak frequency

% Plot settings
iciRange = []; % min/max ici in ms
dbRange = [];  % min/max db for plots of pp and rms
durRange = []; % min/max duration in us
frRange = [fLow, fHi];   % min/max frequency for plots of peak and center freq
frdbwRange = [fLow, fHi]; % min/max frequency for plots of 3/10 db bw
durstep = 1; % step range for number bins in histogram 





