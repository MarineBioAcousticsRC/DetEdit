function writeSFilesTimes(cDir,saveDir)

secInDay = 60*60*24;

% get folder names
folders = dir(cDir) ;
folders = ({folders.name})' ;
folders = folders(3:end);

% get deployment names
siteDiskList = cell(length(folders),1);
for s = 1: length(folders)
    a = cell2mat(strfind(folders(s),'_disk'));
    siteDiskList{s,1} = folders{s}(1:a-1);
end

siteDisk = unique(siteDiskList);

for i = 1:length(siteDisk)
    % loop over files from the same site (e.g., MC01)
    disp(['loading times from: ', siteDisk{i}]);
    
    index = strfind(siteDiskList, siteDisk{i});
    siteDiskIdx = find(not(cellfun('isempty', index)));
    
    % return detection times from the same site
    times = [];
    allfolders = [];
    for j = 1:length(siteDiskIdx)
        f = siteDiskIdx(j);
        fold = fullfile(cDir,folders{f});
        disp(['   ', folders{f}]);
        SearchFileMask = {'*.s'};
        SearchPathMask = {fold};
        SearchRecursiv = 1;
        
        [PathFileList, FileList, ~] = ...
            utFindFiles(SearchFileMask, SearchPathMask, SearchRecursiv);
        
        for idx = 1:length(PathFileList)
            [sNoise, eNoise, Labels] = ioReadLabelFile(PathFileList{idx});
            if ~isempty(sNoise)
                detShip = strmatch('ship', Labels);
                Starts = sNoise(detShip);
                Stops = eNoise(detShip);
                
                undscr = strfind(FileList{idx},'_');
                t = FileList{idx}(undscr(end-1)+1:end-2);
                timeNum = datenum(t,'yymmdd_HHMMSS');
                StartTimes = timeNum + Starts/secInDay;
                StopTimes = timeNum + Stops/secInDay;
                
                times = [times;[StartTimes, StopTimes]];
                allfolders = [allfolders;repmat(folders(siteDiskIdx(j)),length(StartTimes),1)];
            end
        end
    end
    

    % write times in mat file
    FileName = [siteDisk{i},'_',allfolders{1}(end-5:end),'-',allfolders{end}(end-1:end),'.mat'];
    save(fullfile(saveDir,FileName),'times');
end

disp('Done writing ship file times')
