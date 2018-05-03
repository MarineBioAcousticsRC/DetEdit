function [zFD, zID,zMD,fNameList]= buildLabelFiles(matchingFile, sdir)
% initialize variables
zFD = [];
zID = [];
zMD = [];

% Name and build false detection file
ffn = strrep(matchingFile,'TPWS','FD');
fNameList.FD = fullfile(sdir,ffn);
AFD = exist(fNameList.FD,'file');
if (AFD ~= 2) % if it doesn't exist, make it
    zFD(1,1) = 1;
    save(fNameList.FD,'zFD');
    disp('Made new FD file');
end

% Name true detection file
tfn = strrep(matchingFile,'TPWS','TD');
fNameList.TD = fullfile(sdir,tfn);
% NOTE: TD file is made elsewhere because it depends on a later variable


% Name and build ID file
idfn = strrep(matchingFile,'TPWS','ID');
fNameList.ID = fullfile(sdir,idfn);
AID = exist(fNameList.ID,'file');
if (AID ~= 2)% if it doesn't exist, make it
    zID = [];
    save(fNameList.ID,'zID');
    disp('Made new ID file');
end

% Name and build ID file
mdfn = strrep(matchingFile,'TPWS','MD');
fNameList.MD = fullfile(sdir,mdfn);
AMD = exist(fNameList.MD,'file');
if (AMD ~= 2)% if it doesn't exist, make it
    zMD = [];
    save(fNameList.MD,'zMD');
    disp('Made new MD file');
end
