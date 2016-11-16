% mkTPWS_perDir.m

% Script takes output from de_detector.m and put it into a format for use in
% detEdit.
% Output:
% A TPWS.m file containing 4 variables:
%   MTT: An Nx2 vector of detection start and end times, where N is the
%   number of detections
%   MPP: An Nx1 vector of recieved level (RL) amplitudes.
%   MSP: An NxF vector of detection spectra, where F is dictated by the
%   parameters of the fft used to generate the spectra and any
%   normalization preferences.
%   f = An Fx1 frequency vector associated with MSP

clearvars

% Setup variables:
baseDir = 'D:\MPmetadata\MP01'; % directory containing de_detector output
outDir = 'D:\MPmetadata\MP01_TPWS'; % directory where you want to save your TTPP file
siteName = 'GofMX_MP01'; % site name, used to name the output file
ppThresh = 120; % minimum RL in dBpp. If detections have RL below this
% threshold, they will be excluded from the output file. Useful if you have
% an unmanageable number of detections.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the output file exists, if not, make it
if ~exist(outDir,'dir')
    fprintf('Creating output directory %s\n',outDir)
    mkdir(outDir)
end
letterCode = 97:122;
dirSet = dir(fullfile(baseDir,[siteName,'*']));
for itr0 = 1:length(dirSet)
    if dirSet(itr0).isdir &&~strcmp(dirSet(itr0).name,'.')&&...
            ~strcmp(dirSet(itr0).name,'..')
        
        letterFlag = 0; % flag for knowing if a letter should be appended to disk name
        inDir = fullfile(baseDir,dirSet(itr0).name);
        fileSet = what(inDir);
        lfs = length(fileSet.mat);
        clickTimesVec = [];
        ppSignalVec = [];
        specClickTfVec = [];
        tsVecStore = [];
        subTP = 1;

        for itr2 = 1:lfs
            thisFile = fileSet.mat(itr2);
            
            load(char(fullfile(inDir,thisFile)),'-mat','clickTimes','hdr',...
                'ppSignal','specClickTfHR','yFiltBuff','fHR','durClick')
            specClickTf = specClickTfHR;
            if ~isempty(clickTimes)
                keepers = find(ppSignal >= ppThresh);
                
                ppSignal = ppSignal(keepers);
                clickTimes = clickTimes(keepers,:);
                
                [~,keepers2] = unique(clickTimes(:,1));
                
                clickTimes = clickTimes(keepers2,:);
                ppSignal = ppSignal(keepers2);
                
                fileStart = datenum(hdr.start.dvec);
                posDnum = (clickTimes(:,1)/(60*60*24)) + fileStart +...
                    datenum([2000,0,0,0,0,0]);
                clickTimesVec = [clickTimesVec; posDnum];
                ppSignalVec = [ppSignalVec; ppSignal];
                tsWin = 200;
                tsVec = zeros(length(keepers),tsWin);
                f = fHR;
                for iTS = 1:length(keepers(keepers2))
                    thisClick = yFiltBuff{keepers(keepers2(iTS))};
                    [~,maxIdx] = max(thisClick);
                    % want to align clicks by max cycle
                    
                    dTs = (tsWin/2) - maxIdx; % which is bigger, the time series or the window?
                    dTe =  (tsWin/2)- (length(thisClick)-maxIdx); % is the length after the peak bigger than the window?
                    if dTs<=0 % if the signal starts more than N samples ahead of the peak
                        % the start position in the TS vector has to be 1
                        posStart = 1;
                        % the start of the click ts has to be peak - 1/N
                        sigStart = maxIdx - (tsWin/2)+1;
                    else
                        % if it's smaller
                        posStart = dTs+1; % the start has to make up the difference
                        sigStart = 1; % and use the first index of the signal
                    end
                    
                    if dTe<=0 % if it ends after the cut off of N samples
                        posEnd = tsWin; % the end in the TS vector has to be at N
                        sigEnd = maxIdx + (tsWin/2); % and last index is N/2 points after the peak.
                    else % if it ends before the end of the TS vector
                        posEnd = tsWin-dTe; % the end point has to be the
                        % difference between the window length and the click length
                        sigEnd = length(thisClick);
                    end
                    
                    tsVec(iTS,posStart:posEnd) = thisClick(sigStart:sigEnd);
                end
                tsVecStore = [tsVecStore;tsVec];
                
                if iscell(specClickTf)
                    spv = cell2mat(specClickTf');
                    specClickTfVec = [specClickTfVec; spv(:,keepers(keepers2))'];
                else
                    specClickTfVec = [specClickTfVec; specClickTf(keepers(keepers2),:)];
                end
                clickTimes = [];
                hdr = [];
                specClickTf = [];
                ppSignal = [];
            end
            fprintf('Done with file %d of %d \n',itr2,lfs)
            
            if (size(clickTimesVec,1)>= 1800000 && (lfs-itr2>=10))|| itr2 == lfs
                
                MSN = tsVecStore;
                MTT = clickTimesVec;
                MPP = ppSignalVec;
                MSP = specClickTfVec;
                if itr2 == lfs && letterFlag == 0
                    ttppOutName =  [fullfile(outDir,dirSet(itr0).name),'_Delphin_TPWS1','.mat'];
                    fprintf('Done with directory %d of %d \n',itr0,length(dirSet))
                    subTP = 1;
                else
                    
                    ttppOutName = [fullfile(outDir,dirSet(itr0).name),char(letterCode(subTP)),'Delphin_TPWS1','.mat'];
                    subTP = subTP+1;
                    letterFlag = 1;
                end
                save(ttppOutName,'MTT','MPP','MSP','MSN','f','-v7.3')
                
                MTT = [];
                MPP = [];
                MSP = [];
                MSN = [];
                
                clickTimesVec = [];
                ppSignalVec = [];
                specClickTfVec = [];
                tsVecStore = [];
                
            end
        end
    end
end
    
