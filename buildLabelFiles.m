function [zFD, zID,fNameList]= buildLabelFiles(matchingFile, p)

% buildLabelFiles.m

% Takes file name and directory and creates FD and ID files

% Inputs:
%
%   matchingFile - TPWS file name
%
%   p - parameter struct. Optional fields:
%
%       p.labelDir - directory holding the FD/TD/ID label files. Defaults to
%           p.tpwsDir, which is the historical behaviour. Set this to keep
%           labels somewhere other than beside the TPWS files, so several
%           label versions can coexist without copying files around.
%
%       p.labelItr - iteration number used in the FD/TD/ID file names,
%           independent of the TPWS iteration. Defaults to whatever the TPWS
%           file uses. Set this to pair, for example, ID2 with TPWS1: the
%           TPWS file is an immutable input and does not need to be
%           duplicated just to hold a new set of labels.
%
% Outputs:
%
%   zFD - Inicialize variable of detection times to label as false detections
%
%   zID - Inicialize variable of detection times to label as ID detections
%
%   fNameList - A struct with 3 fields indicating the directory path to
%           FD,TD and ID files


zFD = [];
zID = [];

% Where the label files live. Defaults to the TPWS directory.
if isfield(p,'labelDir') && ~isempty(p.labelDir)
    labelDir = p.labelDir;
    if ~exist(labelDir,'dir')
        mkdir(labelDir);
        fprintf('Made label directory: %s\n',labelDir);
    end
else
    labelDir = p.tpwsDir;
end

% Iteration number used in the label file names. Empty means "same as the
% TPWS file", which is the historical behaviour.
if isfield(p,'labelItr') && ~isempty(p.labelItr)
    if isnumeric(p.labelItr)
        labelItr = num2str(p.labelItr);
    else
        labelItr = p.labelItr;
    end
else
    labelItr = '';
end

% Name and build false detection file
ffn = labelFileName(matchingFile,'FD',labelItr);
fNameList.FD = fullfile(labelDir,ffn);
AFD = exist(fNameList.FD,'file');
if (AFD ~= 2) % if it doesn't exist, make it
    zFD(1,1) = 1;
    save(fNameList.FD,'zFD');
    disp('Made new FD file');
end

% Name true detection file
tfn = labelFileName(matchingFile,'TD',labelItr);
fNameList.TD = fullfile(labelDir,tfn);
% NOTE: TD file is made elsewhere because it depends on a later variable


% Name and build ID file
idfn = labelFileName(matchingFile,'ID',labelItr);
fNameList.ID = fullfile(labelDir,idfn);
AID = exist(fNameList.ID,'file');
if (AID ~= 2)% if it doesn't exist, make it
    zID = [];
    mySpID = p.mySpID;

    save(fNameList.ID,'zID','mySpID');
    disp('Made new ID file');
end

if ~strcmp(labelDir,p.tpwsDir) || ~isempty(labelItr)
    fprintf('Labels: %s\n',fNameList.ID);
end

end


function name = labelFileName(matchingFile,token,labelItr)
% Derive a label file name from a TPWS file name.
%
% Replaces the TPWS token, keeping the iteration digits unless labelItr asks
% for different ones. Renumbering is what allows ID2 to sit alongside TPWS1
% instead of requiring a duplicate multi-gigabyte TPWS2.

name = regexprep(matchingFile,'TPWS(\d*)',[token,'$1']);
if ~isempty(labelItr)
    name = regexprep(name,[token,'\d*'],[token,labelItr]);
end

end
