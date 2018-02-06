% itr_mkLTSA.m
% Iterates over a directory of TPWS files and calls mkLTSAsessions for each
% one.

clearvars
clear global
filePrefix = 'SOCAL35S'; % File name to match. 
% File prefix should include deployment, site, (disk is optional). 
% Example: 
% File name 'GofMX_DT01_disk01-08_TPWS2.mat' 
%                    -> filePrefix = 'GofMX_DT01'
% or                 -> filePrefix ='GOM_DT_09' (for files names with GOM)
sp = 'Cuviers'; % your species code
itnum = '2'; % which iteration you are looking for
srate = 200; % sample rate
LTSApath = 'E:\JAH\BW\Detections\SOCAL_Cuviers\SOCAL_S\SOCAL_35_S\SOCAL35S_LTSAs'; % directory containing all LTSAs for this deployment
% LTSA folder should match the site specified in prefix
tpwsPath = 'E:\JAH\BW\Detections\SOCAL_Cuviers\SOCAL_S\SOCAL_35_S\TPWS'; %directory of TPWS files
tfName = 'E:\Harp_TF'; % Directory ...
% with .tf files (directory containing folders with different series ...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum); % iternate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsPath,subfolder));
   end
end

% Find all TPWS files that fit your specifications (does not look in subdirectories)
% Concatenate parts of file name
if isempty(sp)
    detfn = [filePrefix,'.*','TPWS',itnum,'.mat'];
else
    detfn = [filePrefix,'.*',sp,'.*TPWS',itnum,'.mat'];
end
% Get a list of all the files in the start directory
fileList = cellstr(ls(tpwsPath));
% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end

% for each TPWS file found, make LTSA.mat file
for iD = 1:length(fileMatchIdx)
    matchingFile = fileList{fileMatchIdx(iD)};
    detfn = dir(fullfile(tpwsPath,matchingFile));
    
    if exist('tfName','var')
        mkLTSAsessions('filePrefix', filePrefix, 'detfn',detfn.name,...
            'sp', sp, 'lpn', LTSApath, 'sdir', tpwsPath,'srate',srate,'tfName',tfName);
    else
        mkLTSAsessions('filePrefix', filePrefix, 'detfn',detfn.name,...
            'sp', sp, 'lpn', LTSApath, 'sdir', tpwsPath,'srate',srate);
    end
    
end