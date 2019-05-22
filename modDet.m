function modDet(userFunc)
% modDet.m
% modify detections based on analyst input
% originally based on jah/smw GoM beaked whale detector work
%
% 140305 smw
% 140318 jah  141003 jah modified for TPWS files
% 150813 jah for KO
%
% algorithm:
% per site and species
% DT1 = DT0 - zFD
% where
% DT0 = current true detection times
% DT1 = modifed detection times
% from analyst input :
% zFD = false detections times (eg ships, others species)

secInDay = 60*60*24; % convert seconds to days

%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:

% get user input and set up function name
typeInput = exist('userFunc','var');
if typeInput ~= 1
    [userfile,userpath] = uigetfile('*.m',...
        'Select Script with your Data Parameter Settings');
    addpath(userpath) % it adds user folder path to the beggining of the set path
    userFunc = str2func(['@',userfile(1:end-2)]);
end

p = getParams(userFunc,'analysis','modDet');

%% Define subfolder that fit specified iteration
if p.iterationNum > 1
    for id = 2: str2num(p.iterationNum) % iternate id times according to p.iterationNum
        subfolder = ['TPWS',num2str(id)];
        p.tpwsDir = (fullfile(p.tpwsDir,subfolder));
    end
end

%% Check if TPWS file exists (does not look in subdirectories)
% Concatenate parts of file name
if isempty(p.speName)
    detfn = [p.filePrefix,'.*','TPWS',p.iterationNum,'.mat'];
else
    detfn = [p.filePrefix,'.*',p.speName,'.*TPWS',p.iterationNum,'.mat'];
end
% Get a list of all the files in the start directory
fileList = cellstr(ls(p.tpwsDir));
% Find the file name that matches the p.filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end

%% Define output directory
newiterationNum = num2str(str2double(p.iterationNum)+1);
inTPWS = ['TPWS',p.iterationNum];
outTPWS = ['TPWS',newiterationNum];

outDir = fullfile(p.tpwsDir,outTPWS);
if ~isdir(outDir)
    disp(['Make new folder: ',outDir])
    mkdir(outDir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Handle Transfer Function
if (p.tfSelect > 0) || ~isempty(strfind(p.calcParams,'all'))
    [tf,~,~] = getTransfunc(p.filePrefix, p.tfName,p);
else
    tf = 0;
    disp('No TF Applied')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make new TPWS file (TPWS2) for each TPWS file found
for iD = 1:length(fileMatchIdx)
    
    %% Find files and load data
    % Find if TPWS file exist
    matchingFile = fileList{fileMatchIdx(iD)};
    detfn = dir(fullfile(p.tpwsDir,matchingFile));
    
    fNameList.TPWS = fullfile(p.tpwsDir,detfn.name);
    A1 = exist(fNameList.TPWS,'file');
    if A1 ~= 2
        disp(['Error: File Does Not Exist: ',fNameList.TPWS])
        return
    end
    
    % Load detections
    load(fNameList.TPWS)
    
    if isrow(MTT); MTT = MTT'; end
    if isrow(MPP); MPP = MPP'; end
    
    % Apply tf (if defined) and remove low amplitude detections
    MPP = MPP + tf;
    ib = find(MPP >= p.threshRL);
    disp([' Removed too low:',num2str(length(MPP)-length(ib))]);
    MTT = MTT(ib);
    MPP = MPP(ib);
    MSN = MSN(ib,:); 
    MSP = MSP(ib,:); 
    
    disp(['Number of Starting Detections = ',num2str(length(MTT))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Remove false detections and ID detections if specified
    % Load False Detections file
    zFDfn = strrep(detfn.name,inTPWS,['FD',p.iterationNum]);
    load(fullfile(p.tpwsDir,zFDfn))
    
    % Remove False Detections
    [MTT,IA] = setdiff(MTT,zFD'); % setdiff already sorts the data
    MPP = MPP(IA);
    MSN = MSN(IA,:);
    MSP = MSP(IA,:);
    
    % Load ID Detections file
    zIDfn = strrep(detfn.name,inTPWS,['ID',p.iterationNum]);
    load(fullfile(p.tpwsDir,zIDfn))
    
    % Remove ID Detections (if specified)
    if p.excludeID
        [MTT,IB] = setdiff(MTT,zID); % setdiff already sorts the data
        MPP = MPP(IB);
        MSN = MSN(IB,:);
        MSP = MSP(IB,:);
    end
    
    disp(['Number of Final Detections = ',num2str(length(MTT))]);
    
    % Save ID Detections to Output File
    outFileTPWS = strrep(detfn.name,inTPWS,outTPWS);
    outFileID = strrep(outFileTPWS,'TPWS','ID');
    
    disp(['Save ',fullfile(outDir,outFileID)])
    save(fullfile(outDir,outFileID),'zID','-v7.3')
    disp('Done Modifying File')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Check if there is at least one encounter longer than minBout to store file
    if ~ isempty(MTT)
        
        % Calculate bout duration and find encounters longer than minBout
        [~,~,~,bd] = calculate_bouts(MTT,p);
        bdI = find(bd > (p.minBout / secInDay), 1);
        if tf > 0
           MPP = MPP - tf; 
        end
        
        if ~isempty(bdI)
            disp(['Save ',fullfile(outDir,outFileTPWS)])
            if exist('f','var')
                save(fullfile(outDir,outFileTPWS),'f','MTT','MPP','MSN','MSP','-v7.3')
            else
                warning('no frequency vector available')
                save(fullfile(outDir,outFileTPWS),'MTT','MPP','MSN','MSP','-v7.3')
            end
            
            % Calculate parameters and make figure if specified by user in itr_modDet
            if p.calcParams
                paramsHist(MTT,MPP,MSP,outDir,outFileTPWS,p)
            else
                disp('No parameters calculated')
            end
        else
            disp(['No encounter longer than minimum bout (',num2str(p.minBout),') in file: ',outFileTPWS])
        end
    else
        disp(['No true detections in file: ',outFileTPWS])
    end
     
end
