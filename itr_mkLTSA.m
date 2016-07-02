% itr_mkLTSA.m
% Iterates over a directory of TPWS files and calls mkLTSAsessions for each
% one.

clearvars

% setup info - modify to fit your site, species and path
stn = 'DT'; % sitename
dpn = '03'; % deployment code
sp = 'De'; % your species code
itnum = 1; % which iteration you are looking for
LTSApath = 'G:\ltsa\DT\DT03'; % directory containing all LTSAs for this deployment
tpwsPath = 'F:\detStore\DTmetadata\DT03'; %directory of TPWS files



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Find all TPWS files that fit your specifications (does not look in subdirectories)
diskList = dir(fullfile(tpwsPath,[stn,dpn,'_*TPWS',num2str(itnum),'.mat']));

% for each TPWS file found, make LTSA.mat file
for iD = 1:length(diskList);
    mkLTSAsessions('stn', stn, 'dpn', dpn, 'itnum', num2str(itnum),...
       'sp', sp, 'lpn', LTSApath, 'sdir', tpwsPath)
end