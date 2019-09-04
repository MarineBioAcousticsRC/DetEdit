% detEdit_settings
% inputs can be defined here instead of through gui windows

% You can make different versions of this, with different names
% for different species or sites.
%
% eg. detEdit_settings_kogia.m
%
% Then change line 22 of detEdit to match the settings file you want to
% use.

filePrefix = 'HAWAII07_disk04a'; % File name to match.
% File prefix should include deployment, site, (disk is optional).
% Example:
% File name 'GofMX_DT01_disk01-08_TPWS2.mat'
%                    -> filePrefix = 'GofMX_DT01_disk01-08'
iterationNum = '1'; % iteration
sampleRate = 200; % sample rate
sp = ''; % species code (can be: 'Ko' or 'k' (kogia);
% 'Zc' or 'z' (Cuvier's),'Me' or 'm' (Gervais'), 'Md' (Blainville's), BWG,...
% 'De' (Dolphin), 'Po' (porpoise), 'MFA', 'whs' (whistles), 'Dl' (beluga)
c4fd = 100; % Interval to check for false detections
tpwsDir = 'F:\Hawaii_K\TPWS'; %Directory with TPWS files
outdir = 'F:\Hawaii_K\TPWS\detEdit_image\HI07'; %out directory where you want to save images of bouts, if desired
tfName = [];%'E:\TF_files'; % Directory ...
% with .tf files (directory containing folders with different series ...
% (e.g. 300_series,400_series)

% Colors to use for classification
colorTab = [255, 153, 200; ... % type 1 pink
    122,  15, 227; ... % type 2 dark purple
    174, 235, 255; ... % type 3 pale-blue
    0, 255, 255; ... % type 4 cyan
    255, 177, 100; ... % type 5 peach
    255,   0, 255; ... % type 6 magenta
    20,  43, 140; ... % type 7 dark blue
    218, 179, 255; ... % type 2 purple
    179, 200, 255; ... % type 3 light-blue
    221, 125,   0]./255; % type 10  orange
colorTab = round(colorTab.*100)/100;

%settings for screengrab of all plots
%imgDim = [0,0,1450,850];%make sure to check the first couple plots 
%you do and fiddle with this as necessary

%settings for if you want a template spectra included
plotTemplate = 1; %set to 0 if no template desired


if plotTemplate
    load('GmPc_template.mat'); %load in file of templates
    template_color = [0.4 0.4 0.4]; %default is dark gray
    template_color2 = [0.7 0.7 0.7]; %comment this in if using multiple templates. current color is light gray.
    %which template do you want to use?- comment it in
    template_spectra = Gm_specClickTf;
    template_spectra2 = Pc_specClickTf;
end


%% Settings preferences to override defaults
% Comment these in as needed to override detEdit defaults

paramsUser.ltsaLims = [0,100]; % min and max ylimits in kHz for ltsa plot
paramsUser.ltsaMax = 6; % ltsa maximum duration per session
% paramsUser.tfSelect = 0; % freq used for transfer function, leave at 0 if no adjustment
% paramsUser.specChar = 'Unk';  %Simone abbreviation for species
paramsUser.speName = '';  % Species code used in file names
% paramsUser.dtHi = .3; % max yaxis value for IPI display in sec
paramsUser.fLow = 5; % Minimum frequency of interest in kHz
paramsUser.fHi = 60;%  Maximum frequency of interest in kHz

% paramsUser.threshRL = 0; % minimum RL threshold in dB peak-to-peak
% paramsUser.threshRMS = 126; % RMS threshold cutoff
% paramsUser.threshHiFreq = 30; % high freq cutoff for clicks
% paramsUser.ltsaContrast = 116; % ltsa contrast
paramsUser.ltsaBright = 30; % ltsa brightness
paramsUser.ltsaLims = [0,60]; % max and min of LTSA plot
paramsUser.rlLow = 120; % PP plot window low limit
paramsUser.rlHi = 200; % PP plot window high limit
% paramsUser.dfManual = []; % LTSA step size in 10 [Hz] bins
% paramsUser.p1Low = thresRL - 5;
% paramsUser.p1Hi = 170;
paramsUser.minBout = -1; % minimum bout duration in seconds

%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
gth = .5;    % gap time in hrs between sessions
minNdet = 1; % minimum number of detections per session. Sessions with fewer than this will be skipped
maxDetLoad = 1e5; % [] - read all or 4e5 - the number of detections above
% which you want to read from disk instead of loading all spectra and
% timeseries into memory this is for large files (e.g. dolphin click detections)
% if maxDetLoad exist, plotaxes can be defined to keep the format of plot
% 51 and 53
%plotaxes.minRMS = 60;
%plotaxes.maxRMS = 130;
%plotaxes.minPP = 120;
%plotaxes.maxPP = 180;
