% dolphinSettings

% Example script to define directories and parameter preference for the
% interface. You can make different versions of this, with different names
% for different species or sites.

% Define input/output locations.REQUIRED
filePrefix = 'Hawaii17K_newTF_disk01b_Delphin'; % TPWS file name to match. 
% Optional, replace file prefix to a more generic name to specify settings 
% for mkLTSAsessions or modDet, it will run in multiple files.
% Example: GofMX_DT03 (will run modDet to all files matching the generic name) 
iterationNum = '1'; % iteration number
sampleRate = 200; % replace with your sample rate
sp = []; % species code 
% Example:  '' (Unknown), 'De' (Dolphin), 'Pm' (sperm whale)
% (See comments in initSpParams.m for more species codes)
tpwsDir = 'G:\reProcessed\Hawaii\Hawaii17K\TPWS_newTF'; % identify folder containing TPWS files

% Specific input/output locations (comment them if not in use)
% tfName = 'E:\MyTFfolder'; % identify folder containing transfer function 
% files (.tf). Required if spectra has not been calculated peak to received levels 

ltsaDir = 'J:\USWTR02B\LTSAs'; % identify folder containing 
% ltsa files (.ltsa). REQUIRED to run mkLTSAsessions.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setting preferences to override defaults parameters for the interface
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m and initSpParams.m)

%% Bout preferences
% paramsUser.threshRL = 115; % minimum RL threshold in dB peak-to-peak
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
paramsUser.ltsaContrast = 116; % ltsa contrast
paramsUser.ltsaBright = 55; % ltsa brightness
% paramsUser.ltsaLims = [0,sampleRate/2]; % max and min of LTSA plot
paramsUser.ltsaMax = 1; % ltsa maximum duration per session
paramsUser.dtHi = 0.5; % max yaxis value for ICI display in sec
% paramsUser.minDur = []; % minimum window duration (if specified in minutes)

%% Panel Frequency spectra
paramsUser.fLow = 10; % Minimum frequency of interest
paramsUser.fHi = 90; % Maximum frequency of interest

%% Panel RL rms vs. RL pp | Peak freq.
% paramsUser.slope = 0.7; % slope for shifting data vertically
% paramsUser.threshRMS = 0; % default for < command, RMS threshold cutoff
% paramsUser.threshPP = 0; % default for : command, PP threshold cutoff
% paramsUser.threshHiFreq = 40; % default for ^ command, high freq cutoff for clicks

paramsUser.rmsLow = -40;
paramsUser.rmsHi = 40;

%%% Setting preferences to override modDet.m defaults
% (Uncomment these in as needed to override detEdit defaults defined at
% initDefaultParams.m)

% paramsUser.excludeID = 1; % yes - 1 | no - 0. Exclude ID times from MTT files 
% paramsUser.calcParams = 1;% yes - 1 | no - 0. Calculate Parameters peak-to-peak, 
% inter-detection-interval and peak frequency

%% Other settings
paramsUser.nTestBins = 25; %number of bins for fpfn tests 
paramsUser.minClicks = 10; %minimum number of clicks to include a bin for testing, can set to 0 if not needed 
% Colors to use for classification - ID signal types
colorTab = [204, 204, 255; ... % Blainville's, lilac
            255, 128, 0; ... % Cuvier's, orange
            102, 255, 178; ... % FKW, seafoam
            0, 153, 0; ... % LF1, crayola green
            0, 128, 255; ... % SFPW1, medium blue
            51, 255, 51; ... % SFPW2, bright green
            25,25,112; ... % Sten1/2, dark blue 
%                        102, 178, 255; ... % Sten1/2, periwinkle
            0, 255, 255; ... % Sten3, cyan
                        255, 0, 255; ... % bott/MHW, magenta
%             255, 204, 153; ... % bott/MHW, tan
            245, 194, 66;... % kogia, mustard yellow
            249, 177, 211]./255; % noise, light pink
        %             255, 0, 0;... % HFA, red
%             153, 51, 255; ... % Kogia, purple
%             255, 0, 0;... % MFA, red
%             255, 0, 0;... % MultiFreq_Sonar, red
%             204, 255, 153; ... % Rissos, lime green
%             255, 0, 0;... % Snapping Shrimp, red
%             102, 51, 0; ... % Sperm whale, brown
%             255, 0, 0;... % , red
paramsUser.colorTab = round(colorTab.*100)/100;      

paramsUser.mySpID = struct('Name',{'Blainvilles','Cuviers','FKW','LF1','SFPW1','SFPW2','Sten1_2',...
    'Sten3','bott_MHW','kogia','noise_sperm'},'zID_Label',{1,2,3,4,5,6,7,8,9,10,11});

