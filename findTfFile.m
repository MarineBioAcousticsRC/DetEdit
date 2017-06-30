function PathFileList = findTfFile(indir,stationDeploy)


% select which tf corresponds to station/deployment
switch stationDeploy
    % site DT
    case {'DT01','DT02','DT03'} 
        tfnum = 589;
        serie = '500_series';
    case {'DT04','DT05','DT06','DT07','DT08'} 
        tfnum = 638;
        serie = '600_series';
    case 'DT09'
        tfnum = 715;
        serie = '700_series';
    
    % site MC
    case {'MC01','MC02','MC03','MC04','MC05'} 
        tfnum = 585;
        serie = '500_series';
    case 'MC06'
        tfnum = 651;
        serie = '600_series';
    case 'MC07'
        tfnum = 693;
        serie = '600_series';
    case 'MC09'
        tfnum = 729;
        serie = '700_series';
    case 'MC10'
        tfnum = 718;
        serie = '700_series';
    case 'MC11'
        tfnum = 618;
        serie = '600_series';
        
    % site GC
    case {'GC01','GC02','GC03','GC04'} 
        tfnum = 601; 
        serie = '600_series';
    case 'GC05'
        tfnum = 656; 
        serie = '600_series';
    case {'GC06','GC07'}
        tfnum = 694;
        serie = '600_series';
    case 'GC08'
        tfnum = 719;
        serie = '700_series';
    case 'GC09'
        tfnum = 718; 
        serie = '700_series';
        
    otherwise
        tfnum = [];
        disp('Transfer function folder not found or site matches')
end

pathSeries = fullfile(indir,serie);
folders = dir(pathSeries);
folders = ({folders.name})' ;
folders = folders(3:end);

foldNums = cell2mat(strtok(folders,'_')); % get part of number
tfnumFold = find(str2num(foldNums) == tfnum,1,'last'); 

pathTfFile = fullfile(pathSeries,folders{tfnumFold});

SearchFileMask = {'*.tf'};
SearchPathMask = {pathTfFile};
SearchRecursiv = 1;

[PathFileList, ~, ~] = ...
    utFindFiles(SearchFileMask, SearchPathMask, SearchRecursiv);

PathFileList = cell2mat(PathFileList);













