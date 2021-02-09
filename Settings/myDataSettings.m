% yourDataSettings

% Example script to define directories and parameter preference for the
% interface. You can make different versions of this, with different names
% for different species or sites.

% Define input/output locations.REQUIRED
filePrefix = 'Filename'; % TPWS file name to match. 
% Optional, replace file prefix to a more generic name to specify settings 
% for mkLTSAsessions or modDet, it will run in multiple files.
% Example: GofMX_DT03 (will run modDet to all files matching the generic name) 
iterationNum = '1'; % iteration number
sampleRate = 200; % replace with your sample rate
sp = ''; % species code 
% Example:  '' (Unknown), 'De' (Dolphin), 'Pm' (sperm whale)
% (See comments in initSpParams.m for more species codes)
tpwsDir = 'E:\MyTPWSfolder'; % identify folder containing TPWS files

% Specific input/output locations (comment them if not in use)
% tfName = 'E:\MyTFfolder'; % identify folder containing transfer function 
% files (.tf). Required if spectra has not been calculated peak to received levels 

% REQUIRED to run mkLTSAsessions.m:
ltsaDir = 'E:\MyLTSAfolder'; % identify folder containing ltsa files (.ltsa)

% REQUIRED to run summaryParams.m:
effortTimes = 'E:\Effort.xls'; % specify excel file with effort times
referenceTime = '2010-04-01'; %reference time format 'yyyy-MM-dd'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting preferences to override defaults parameters for the interface
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m and initSpParams.m)

%% Bout preferences
% paramsUser.threshRL = 0; % minimum RL threshold in dB peak-to-peak
% paramsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
% paramsUser.minBout = 0;% minimum bout duration in seconds
% paramsUser.gth = .5;    % gap time in hrs between sessions
% paramsUser.binDur = 5; % bin duration in minutes
% paramsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% paramsUser.specploton = 1; %(1 = yes | 0 = no ) spectra parameters plot
% paramsUser.minNdet = 1; % minimum detections per session.
% paramsUser.maxDetLoad = 4e5; % ([] - all) max detections to read from disk ( recommended for large files)
% paramsUser.c4fd = 1; % Detections step size to estimate false detection rate

%% Panel LTSA and time series
% paramsUser.rlLow = 110; % PP plot window low limit
% paramsUser.rlHi = 170; % PP plot window high limit
% paramsUser.ltsaContrast = 250; % ltsa contrast
% paramsUser.ltsaBright = 100; % ltsa brightness
% paramsUser.ltsaLims = [0,sampleRate/2]; % max and min of LTSA plot
% paramsUser.ltsaMax = 6; % ltsa maximum duration per session
% paramsUser.dtHi = .5; % max yaxis value for ICI display in sec
% paramsUser.minDur = []; % minimum window duration (if specified in minutes)

%% Panel Frequency spectra
% paramsUser.fLow = 0; % Minimum frequency of interest
% paramsUser.fHi = sampleRate/2; % Maximum frequency of interest

%% Panel RL rms vs. RL pp | Peak freq.
% paramsUser.slope = 1; % slope for shifting data vertically
% paramsUser.threshRMS = 0; % default for < command, RMS threshold cutoff
% paramsUser.threshPP = 0; % default for : command, PP threshold cutoff
% paramsUser.threshHiFreq = 0; % default for ^ command, high freq cutoff for clicks

%%% Setting preferences to override modDet.m defaults
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m)

% paramsUser.excludeID = 0; % yes - 1 | no - 0. Exclude ID times from MTT files 
% paramsUser.calcParams = 0;% yes - 1 | no - 0. Calculate Parameters peak-to-peak, 
% inter-detection-interval and peak frequency

%% ID labels legend
% Assign an ID to colors
% paramsUser.mySpID = {'Label1','Label2','Label3', 'Label4', 'Label5', 'Label6', 'Label7', 'Label8', 'Label9', 'Label10'};