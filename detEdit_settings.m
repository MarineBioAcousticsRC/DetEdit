% detEdit_settings
% inputs can be defined here instead of through gui windows

% You can make different versions of this, with different names
% for different species or sites.
%
% eg. detEdit_settings_kogia.m
% 
% Then change line 18 of detEdit to match the settings file you want to
% use.

% TODO: Ideally, species-specific preferences could also be moved into here

stn = 'JAX11D_'; % site name
dpn = 'disk01'; % deployment number
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'De'; % species code
c4fd = 100; %Interval to check for false detections
sdir = 'G:\JAX11D_TPWS_testingDetEdit'; %Directory with TPWS files
% tfName = 'E:\Code\TF_files\TF Files\HARP\400_series\492_090309\492_090309_invSensit.tf';


%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
gth = .5;    % gap time in hrs between sessions
minNdet = 1; % minimum number of detections per bin
maxDetLoad = 4e5; % the number of detections above which you want to 
% read from disk instead of loading all spectra and timeseries into memory
% this is for large files (e.g. dolphin click detections)


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


% Settings preferences to override defaults:
spParamsUser.ltsaLims = [0,100]; % min and max ylimits in kHz for ltsa plot
spParamsUser.ltsaMax = 3; % length of ltsa window



