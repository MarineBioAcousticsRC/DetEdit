% detEdit_settings
% inputs can be defined here instead of through gui windows
% TODO: Ideally, species-specific preferences could also be moved into here

stn = 'SOCAL'; % site name
dpn = '34E'; % deployment number
itnum = '1'; % iteration
srate = 200; % sample rate
sp = 'Zc'; % species code
c4fd = 1000; %Interval to check for false detections
sdir = 'G:\'; %Directory with TPWS files
tfName = 'E:\Code\TF_files\TF Files\HARP\400_series\492_090309\492_090309_invSensit.tf';



%%%%%%%% other preferences - modify with care %%%%%%%%%%
specploton = 1; %1 = yes spec plot 0 = no spec plot
ltsamax = 6; % length of ltsa window
gth = .5;    % gap time in hrs between sessions
maxDetLoad = 1e5; % the number of detections above which you want to 
% read from disk instead of loading all spectra and timeseries into memory
% this is for large files (e.g. dolphin click detections)
