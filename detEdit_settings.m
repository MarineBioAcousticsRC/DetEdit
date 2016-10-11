% detEdit_settings
% inputs can be defined here instead of through gui windows

% You can make different versions of this, with different names
% for different species or sites.
%
% eg. detEdit_settings_kogia.m
% 
% Then change line 22 of detEdit to match the settings file you want to
% use.

stn = 'WAT_HZ_01_'; % site name
dpn = 'disk14'; % deployment number
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'De'; % species code (can be: 'Ko' or 'k' (kogia);
% 'Zc' or 'z' (Cuvier's),'Me' or 'm' (Gervais'), 'Md' (Blainville's), BWG,...
% 'De' (Dolphin), 'Po' (porpoise), 'MFA', 'whs' (whistles), 'Dl' (beluga)
c4fd = 100; % Interval to check for false detections
sdir = 'I:\WAT_HZmetadata\TPWS\ready for detEdit'; %Directory with TPWS files
% tfName = 'E:\Code\TF_files\Recalculated\tf files\707_130408\707_130408_invSensit.tf';

% Colors to use for classification
colorTab = [255, 153, 200; ... % type 1 pink
    218, 179, 255; ... % type 2 purple
    179, 200, 255; ... % type 3 light-blue
    174, 235, 255; ... % type 4 pale-blue
    0, 255, 255; ... % type 5 cyan
    0, 255,   0; ... % type 6 bright green
    255, 177, 100; ... % type 8 peach
    255,   0, 255; ... % type 9 magenta
    122,  15, 227; ... % type 10 purple
    20,  43, 140; ... % type 11 dark blue
    221, 125,   0]./255; % orange
colorTab = round(colorTab.*100)/100;


%% Settings preferences to override defaults
% Comment these in as needed to override detEdit defaults

spParamsUser.ltsaLims = [0,100]; % min and max ylimits in kHz for ltsa plot
spParamsUser.ltsaMax = 6; % ltsa maximum duration per session
spParamsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
% spParamsUser.specChar = 'Unk';  %Simone abbreviation for species
% spParamsUser.speName = 'Unknown';  % Species code used in file names 
spParamsUser.dtHi = .4; % max yaxis value for IPI display in sec
% spParamsUser.fLow = 5; % boundary for spectrum plot
spParamsUser.threshRL = 0; % minimum RL threshold in dB peak-to-peak
spParamsUser.ltsaContrast = 150; % ltsa contrast
spParamsUser.ltsaBright = 55; % ltsa brightness
% spParamsUser.ltsaLims = [0,100]; % max and min of LTSA plot
% spParamsUser.rlLow = 110; % PP plot window low limit
% spParamsUser.rlHi = 170; % PP plot window high limit
% spParamsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% spParamsUser.p1Low = thresRL - 5;
% spParamsUser.p1Hi = 170;


%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
gth = .5;    % gap time in hrs between sessions
minNdet = 1; % minimum number of detections per session. Sessions with fewer than this will be skipped
maxDetLoad = 4e5; % the number of detections above which you want to 
% read from disk instead of loading all spectra and timeseries into memory
% this is for large files (e.g. dolphin click detections)

