% dolphinSettings

% Example script to define directories and parameter preference for the
% interface. You can make different versions of this, with different names
% for different species or sites.

% Define input/output locations.REQUIRED
filePrefix = 'WAT_WC_03_disk02f'; % TPWS file name to match. 
% Optional, replace file prefix to a more generic name to specify settings 
% for mkLTSAsessions or modDet, it will run in multiple files.
% Example: GofMX_DT03 (will run modDet to all files matching the generic name) 
iterationNum = '1'; % iteration number
sampleRate = 200; % replace with your sample rate
sp = []; % species code 
% Example:  '' (Unknown), 'De' (Dolphin), 'Pm' (sperm whale)
% (See comments in initSpParams.m for more species codes)
tpwsDir = 'J:\WAT_WC_03\TPWS'; % identify folder containing TPWS files

% Specific input/output locations (comment them if not in use)
% tfName = 'E:\MyTFfolder'; % identify folder containing transfer function 
% files (.tf). Required if spectra has not been calculated peak to received levels 

ltsaDir = 'J:\WAT_WC_03\LTSAs'; % identify folder containing 
% ltsa files (.ltsa). REQUIRED to run mkLTSAsessions.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting preferences to override defaults parameters for the interface
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m and initSpParams.m)

%% Bout preferences
paramsUser.threshRL = 118; % minimum RL threshold in dB peak-to-peak
% paramsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
paramsUser.minBout = 0;% minimum bout duration in seconds
% paramsUser.gth = .5;    % gap time in hrs between sessions
% paramsUser.binDur = 5; % bin duration in minutes
% paramsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% paramsUser.specploton = 1; %(1 = yes | 0 = no ) spectra parameters plot
% paramsUser.minNdet = 1; % minimum detections per session.
% paramsUser.maxDetLoad = 4e5; % ([] - all) max detections to read from disk ( recommended for large files)
paramsUser.c4fd = 1000; % Min # detections to test to estimate false detection rate
paramsUser.nTestBins = 100;

%% Panel LTSA and time series
paramsUser.rlLow = 115; % PP plot window low limit
paramsUser.rlHi = 165; % PP plot window high limit
paramsUser.ltsaContrast = 118; % ltsa contrast
paramsUser.ltsaBright = 53; % ltsa brightness
% paramsUser.ltsaLims = [0,sampleRate/2]; % max and min of LTSA plot
paramsUser.ltsaMax = 6; % ltsa maximum duration per session
paramsUser.dtHi = 1; % max yaxis value for ICI display in sec
% paramsUser.minDur = []; % minimum window duration (if specified in minutes)

%% Panel Frequency spectra
paramsUser.fLow = 1; % Minimum frequency of interest
% paramsUser.fHi = sampleRate/2; % Maximum frequency of interest

%% Panel RL rms vs. RL pp | Peak freq.
% paramsUser.slope = 0.7; % slope for shifting data vertically
% paramsUser.threshRMS = 0; % default for < command, RMS threshold cutoff
% paramsUser.threshPP = 0; % default for : command, PP threshold cutoff
% paramsUser.threshHiFreq = 40; % default for ^ command, high freq cutoff for clicks

%%% Setting preferences to override modDet.m defaults
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m)

% paramsUser.excludeID = 1; % yes - 1 | no - 0. Exclude ID times from MTT files 
% paramsUser.calcParams = 1;% yes - 1 | no - 0. Calculate Parameters peak-to-peak, 
% inter-detection-interval and peak frequency

%% Other settings

% Colors to use for classification - ID signal types
paramsUser.colorTab = [204, 204, 255; ... % Blainville's, lilac
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
paramsUser.colorTab = round(paramsUser.colorTab.*100)/100;      

paramsUser.mySpID = {'Blainvilles','Boats','CT11','CT2+CT9','CT3+CT7','CT4/6+CT10',...
    'CT5','CT8','Cuviers','Gervais','GoM_Gervais','HFA','Kogia','MFA',...
    'MultiFreq_Sonar','Rissos','Snapping_Shrimp','Sowerbys','Sperm_Whale',...
    'Trues'};