% detEdit_settings
% inputs can be defined here instead of through gui windows

% You can make different versions of this, with different names
% for different species or sites.
%
% eg. detEdit_settings_kogia.m
% 
% Then change line 22 of detEdit to match the settings file you want to
% use.

filePrefix = 'GofMX_GC06_disk12'; % File name to match. 
% File prefix should include deployment, site, (disk is optional). 
% Example: 
% File name 'GofMX_DT01_disk01-08_TPWS2.mat' 
%                    -> filePrefix = 'GofMX_DT01_disk01-08'
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'Pm'; % species code (can be: 'Ko' or 'k' (kogia);
% 'Zc' or 'z' (Cuvier's),'Me' or 'm' (Gervais'), 'Md' (Blainville's), BWG,...
% 'De' (Dolphin), 'Po' (porpoise), 'MFA', 'whs' (whistles), 'Dl' (beluga)
c4fd = 3000; % Interval to check for false detections
sdir = 'E:\TPWS'; %Directory with TPWS files
%tfName = 'E:\TF_files'; % Directory ...
% with .tf files (directory containing folders with different series ...
% (e.g. 300_series,400_series)

% Colors to use for classification
colorTab = [255, 153, 200; ... % type 1 pink
    218, 179, 255; ... % type 2 purple
    179, 200, 255; ... % type 3 light-blue
    174, 235, 255; ... % type 4 pale-blue
    0, 255, 255; ... % type 5 cyan
    255, 177, 100; ... % type 6 peach
    255,   0, 255; ... % type 7 magenta
    122,  15, 227; ... % type 8 purple
    20,  43, 140; ... % type 9 dark blue
    221, 125,   0]./255; % type 10  orange
colorTab = round(colorTab.*100)/100;


%% Settings preferences to override defaults
% Comment these in as needed to override detEdit defaults

spParamsUser.ltsaLims = [0,100]; % min and max ylimits in kHz for ltsa plot
spParamsUser.ltsaMax = 6; % ltsa maximum duration per session
spParamsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
% spParamsUser.specChar = 'Unk';  %Simone abbreviation for species
% spParamsUser.speName = 'Pm';  % Species code used in file names 
% spParamsUser.dtHi = .3; % max yaxis value for IPI display in sec
% spParamsUser.fLow = 5; % Minimum frequency of interest in kHz
% spParamsUser.fHi = 10;%  Maximum frequency of interest in kHz

% spParamsUser.threshRL = 120; % minimum RL threshold in dB peak-to-peak
spParamsUser.threshRMS = 0; % RMS threshold cutoff
spParamsUser.threshHiFreq = 30; % high freq cutoff for clicks
spParamsUser.threshPP = 0; % high freq cutoff for clicks
%spParamsUser.slope = 0.8;
% % spParamsUser.ltsaContrast = 116; % ltsa contrast
% spParamsUser.ltsaBright = 55; % ltsa brightness
% % spParamsUser.ltsaLims = [0,100]; % max and min of LTSA plot
% spParamsUser.rlLow = 115; % PP plot window low limit
% spParamsUser.rlHi = 165; % PP plot window high limit
% spParamsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% spParamsUser.p1Low = thresRL - 5;
% spParamsUser.p1Hi = 170;
% spParamsUser.minBout = 75; % minimum bout duration in seconds

%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
gth = .5;  % gap time in hrs between sessions
minNdet = 5; % minimum number of detections per session. Sessions with fewer than this will be skipped
maxDetLoad = 4e5; % [] - read all or 4e5 - the number of detections above 
% which you want to read from disk instead of loading all spectra and 
% timeseries into memory this is for large files (e.g. dolphin click detections)
% if maxDetLoad exist, plotaxes can be defined to keep the format of plot
% 51 and 53
plotaxes.minRMS = 60;
plotaxes.maxRMS = 130;
plotaxes.minPP = 120;
plotaxes.maxPP = 180;
