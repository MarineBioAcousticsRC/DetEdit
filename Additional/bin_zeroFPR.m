%bin_zeroFPR.m

% Use this code if you know certainly that you do not have any false bin in
% your dataset. This will set the False Bin Rate to zero

% Define input/output locations
inDir = 'E:\MyFolder'; % identify folder containing TD files
myFileFlag = '.*TD1.mat'; % include a string to match for identifying the files to be processed.


% Get a list of all the files in the start directory
fileList = cellstr(ls(inDir));

% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,myFileFlag))>0);
TDfiles = fileList(fileMatchIdx);

% Load file and replace false bin value to zero
for i = 1:length(TDfiles)
   load(fullfile(inDir,TDfiles{i}))
   
   check = find(zTD(:,1) ~= -1);
   zTD(check,3:4) = 0;
    
   save(fullfile(inDir,TDfiles{i}),'zTD');
   sprintf([TDfiles{i},' saved'])
end