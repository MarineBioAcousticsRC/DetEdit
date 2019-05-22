function [zFD, zID,fNameList]= buildLabelFiles(matchingFile, p)

% buildLabelFiles.m

% Takes file name and directory and creates FD and ID files

% Inputs:
%
%   matchingFile - TPWS file name
%
%   sdir - Directory path containing TPWS files
%
%
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

% Name and build false detection file
ffn = strrep(matchingFile,'TPWS','FD');
fNameList.FD = fullfile(p.tpwsDir,ffn);
AFD = exist(fNameList.FD,'file');
if (AFD ~= 2) % if it doesn't exist, make it
    zFD(1,1) = 1;
    save(fNameList.FD,'zFD');
    disp('Made new FD file');
end

% Name true detection file
tfn = strrep(matchingFile,'TPWS','TD');
fNameList.TD = fullfile(p.tpwsDir,tfn);
% NOTE: TD file is made elsewhere because it depends on a later variable


% Name and build ID file
idfn = strrep(matchingFile,'TPWS','ID');
fNameList.ID = fullfile(p.tpwsDir,idfn);
AID = exist(fNameList.ID,'file');
if (AID ~= 2)% if it doesn't exist, make it
    zID = [];
    mySpID = p.mySpID;
    save(fNameList.ID,'zID','mySpID');
    disp('Made new ID file');  
end
