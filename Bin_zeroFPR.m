% Zero False Bin Rate
% Use this code if you know certainly that you do not have any false bin in
% your dataset. This will set the False Bir Rate to zero

% select directory where TD files are located 
filesDir = 'E:\TPWS\TPWS2\TPWS3';
filePrefix = '.*TD3.mat';


% Get a list of all the files in the start directory
fileList = cellstr(ls(filesDir));
% Find the file name that matches the filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,filePrefix))>0);

TDfiles = fileList(fileMatchIdx);

for i = 1:length(TDfiles)
   load(fullfile(filesDir,TDfiles{i}))
   
   check = find(zTD(:,1) ~= -1);
   zTD(check,3:4) = 0;
    
   save(fullfile(filesDir,TDfiles{i}),'zTD');
   sprintf([TDfiles{i},' saved'])
end