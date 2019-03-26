%initDefaultParams

% Initialize default parameters for the interface, mkLTSAsessions and modDet

if ~exist('sampleRate','var')
    sampleRate = 200;
    disp('Default Sample Rate: 200 kHz. Modify in your settings script')
end
speName = '';
ltsaDir = '';
tfName = '';

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
threshRMS = 0; % default for < command, RMS threshold cutoff
threshPP = 0; % default for : command, PP threshold cutoff
threshHiFreq = 0; % default for ^ command, high freq cutoff for clicks

% Colors to use for classification - ID signal types
colorTab = [191, 191, 0; ... % type 1 green
            191, 0, 191; ... % type 2 purple
            0, 127, 0; ... % type 3 dark-green
            0, 191, 191; ... % type 4 light-blue
            20, 43, 140; ... % type 5 dark-blue
            218, 179, 255; ... % type 6 pale-purple
            255,   214, 0; ... % type 7 yellow
            222,  125, 0; ... % type 8 orange
            255,  153, 199; ... % type 9 pink
            153, 51,   0]./255; % type 10  brown
colorTab = round(colorTab.*100)/100;      


% Parameters for modDet (if specified to compute)

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



