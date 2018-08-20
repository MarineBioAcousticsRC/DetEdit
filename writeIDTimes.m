function writeIDTimes(cDir,saveDir,maxItr)
% Specify directory where to save reduced mat files
sec2Day = 60*60*24;
addms = 1/(1000*sec2Day); % convert 1ms to days
encDur = 300; % 5min apart (300s)

% Get a list of all the files in the start directory
fileList = cellstr(ls(cDir));
% Find the file name that matches the filePrefix
detfn = ['.*ID.*','.mat'];
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
TPWSList = fileList(fileMatchIdx);

siteDiskList = {};
for s = 1: length(TPWSList)
    a = cell2mat(strfind(TPWSList(s),'_disk'));
    siteDiskList{s,1} = TPWSList{s}(1:a-1);
end

siteDisk = unique(siteDiskList);
expression = 'ID(\d+).';
for i = 1:length(siteDisk)
    % loop over files from the same site (e.g., MC01)
    disp(['loading times from: ', siteDisk{i}]);
    
    index = strfind(siteDiskList, siteDisk{i});
    siteDiskIdx = find(not(cellfun('isempty', index)));
    
    % return detection times from the same site
    timesStart = [];
    for j = 1:length(siteDiskIdx)
        f = siteDiskIdx(j);
        file = TPWSList{f};
        fold = fullfile(cDir,file);
        itr = 1;
        cDirTemp = cDir;
        while itr < maxItr+1
            if itr ~= 1
               subfolder = ['TPWS',num2str(itr)];
               cDirTemp = fullfile(cDir,subfolder);
               replace = ['ID',num2str(itr),'.'];
               file = regexprep(file,expression,replace);
               fold = fullfile(cDirTemp,file);
            end
            disp(['   ', file]);
            if exist(fold,'file')
                load(fold);
                if ~isempty(zID)
                    timesStart = [timesStart;zID(:,1)];
                end
                itr = itr+1;
            else
                disp(['   ',file,' does not exist '])
                itr = maxItr+1;
            end
        end
    end
    
    if ~isempty(timesStart)
    timesStart = sort(timesStart);
    diffStarts = diff(timesStart);
    breaks = find(diffStarts>(encDur/sec2Day)); 
    starts = [timesStart(1); timesStart(breaks+1)]; % start 1ms after
    ends = [timesStart(breaks)+datenum(0,0,0,0,0,1/1000); timesStart(end)+datenum(0,0,0,0,0,1/1000)];
    times = [starts,ends];
    % write times in mat file
    FileName = [siteDisk{i},'.mat'];
    save(fullfile(saveDir,FileName),'times');
    end
end

disp('Done')
