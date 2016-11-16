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

stn = 'DT'; % site name
dpn = '02_disk01'; % deployment number
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'De'; % species code
c4fd = 1000; %Interval to check for false detections
sdir = 'G:\DTmetadata\DT02\TPWS'; %Directory with TPWS files
% tfName = 'E:\Code\TF_files\TF Files\HARP\400_series\492_090309\492_090309_invSensit.tf';


%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
ltsamax = 1; % max length of ltsa window
gth = .5;    % gap time in hrs between sessions
minNdet = 1; % minimum number of detections per bin
maxDetLoad = 1e5; % the number of detections above which you want to 
% read from disk instead of loading all spectra and timeseries into memory
% this is for large files (e.g. dolphin click detections)
