% detEdit_settings
% inputs can be defined here instead of through gui windows

% You can make different versions of this, with different names
% for different species or sites.
%
% eg. detEdit_settings_kogia.m
% 
% Then change line 22 of detEdit to match the settings file you want to
% use.

filePrefix = 'GofMX_DT02_disk01'; % File name to match. 
% File prefix should include deployment, site, (disk is optional). 
% Example: 
% File name 'GofMX_DT01_disk01-08_TPWS2.mat' 
%                    -> filePrefix = 'GofMX_DT01_disk01-08'
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'Pm'; % species code (can be: 'Ko' or 'k' (kogia);
% 'Zc' or 'z' (Cuvier's),'Me' or 'm' (Gervais'), 'Md' (Blainville's), BWG,...
% 'De' (Dolphin), 'Po' (porpoise), 'MFA', 'whs' (whistles), 'Dl' (beluga)
sdir = 'C:\TPWS'; %Directory with TPWS files
tfName = ''; %'E:\TF_files'; % Directory ...
% with .tf files (directory containing folders with different series ...
% (e.g. 300_series,400_series)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Settings preferences to override defaults parameters for the interface
% Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams

%% Bout preferences
% paramsUser.threshRL = 0; % minimum RL threshold in dB peak-to-peak
% paramsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
% paramsUser.minBout = 0;% minimum bout duration in seconds
% paramsUser.p.gth = .5;    % gap time in hrs between sessions
% paramsUser.binDur = 5; % bin duration in minutes
% paramsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% paramsUser.specploton = 1; %(1 = yes | 0 = no ) spectra parameters plot
% paramsUser.p.minNdet = 1; % minimum detections per session.
% paramsUser.maxDetLoad = 4e5; % ([] - all) max detections to read from disk ( recommended for large files)
% paramsUser.c4fd = 1; % Detections step size to estimate false detection rate

%% Panel LTSA and time series
% paramsUser.rlLow = 110; % PP plot window low limit
% paramsUser.rlHi = 170; % PP plot window high limit
% paramsUser.ltsaContrast = 250; % ltsa contrast
% paramsUser.ltsaBright = 100; % ltsa brightness
% paramsUser.ltsaLims = [0,srate/2]; % max and min of LTSA plot
% paramsUser.ltsaMax = 6; % ltsa maximum duration per session
% paramsUser.dtHi = .5; % max yaxis value for ICI display in sec
% paramsUser.minDur = []; % minimum window duration (if specified in minutes)

%% Panel Frequency spectra
% paramsUser.fLow = 0; % Minimum frequency of interest
% paramsUser.fHi = srate/2; % Maximum frequency of interest

%% Panel RL rms vs. RL pp | Peak freq.
% paramsUser.slope = 1; % slope for shifting data vertically
% paramsUser.threshRMS = 0; % default for < command, RMS threshold cutoff
% paramsUser.threshPP = 0; % default for : command, PP threshold cutoff
% paramsUser.threshHiFreq = 0; % default for ^ command, high freq cutoff for clicks

%% Colors to use for classification - ID signal types
%{
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
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
