% itr_modDet
% Iterates over a directory of TPWS files and calls modDet for each
% one.

clearvars
clear global
filePrefix = 'GofMX_GC01'; % File name to match. 
% File prefix should include deployment, site, (disk is optional). 
% Example: 
% File name 'GofMX_DT01_disk01-08_TPWS2.mat' 
%                    -> filePrefix = 'GofMX_DT01'
% or                 -> filePrefix ='GOM_DT_09' (for files names with GOM)
sp = 'Pm'; % your species code
itnum = '1'; % which iteration you are looking for
getParams = 'ici&pp'; % Calculate Parameterss: 
%                   -> 'none' do NOT compute parameters
%                   -> 'ici&pp' only to compute peak-to-peak and ici
%                   -> 'all' compute pp, ici, 3/10dbBw, peakFr, F0, rms, dur
srate = 200; % sample rate
gth = .5;  % gap time in hrs between sessions
tpwsPath = 'E:\TPWS'; %directory of TPWS files
%tfName = 'E:\transfer_functions'; % Directory ...
% with .tf files (directory containing folders with different series ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

% for each TPWS file found, make TPWS(itnum+1).mat file
for iD = 1:length(fileMatchIdx);
    matchingFile = fileList{fileMatchIdx(iD)};
    detfn = dir(fullfile(tpwsPath,matchingFile));
    
    if exist('tfName','var')
    modDet('filePrefix', filePrefix, 'detfn',detfn.name,...
       'sp', sp, 'sdir', tpwsPath,'srate',srate,'itnum', itnum,...
       'getParams',getParams,'tfName',tfName)
    else
        modDet('filePrefix', filePrefix, 'detfn',detfn.name,...
            'sp', sp, 'sdir', tpwsPath,'srate',srate,'itnum', itnum,...
            'getParams',getParams)
    end
end

disp('Done processing')