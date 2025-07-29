function [allNewSpec, allNewWaveEnv, allNewICI, allNewID] = make_corrected_bin_clusters(zIDcorrected, zIDorig, TPWSName)
binSize = 5;% bin size in minutes
minClick = 5; % only make clusters for bins with this number of clicks or more
ICIvec = 0:.01:1;
% Find events that were labeled differently in zIDcorrected, excluding
% things with no label.

[C,iCorr,iOrig] = intersect(zIDcorrected(:,1),zIDorig(:,1));


updatedSet = find(zIDcorrected(iCorr,2)~=zIDorig(iOrig,2));
if isempty(updatedSet)
    disp('no corrections to make, skipping ID file')
    return
end
% if there are corrections, we will need the TPWS data
load(TPWSName,'MTT')
TPWSfileObj = matfile(TPWSName);
wMSP = size(TPWSfileObj,'MSP',2);
wMSN = size(TPWSfileObj,'MSN',2);

setToBin = zIDcorrected(updatedSet,:);

timeBins = (floor(min(setToBin(:,1))):(1/24/60)*binSize:ceil(max(setToBin)))';
[Ch,Ih] = histc(setToBin(:,1),timeBins);

%%%timeBins(Ch==0,:) = [];
allNewSpec = [];
allNewWaveEnv = [];
allNewICI = [];
allNewID = [];
for iBin = 1:length(timeBins)
    if Ch(iBin)>minClick
        thisClickSet = setToBin(Ih==iBin,:);
        uSp = unique(thisClickSet(:,2));
        nSp = length(uSp);
        for iSp = 1:nSp
            thisClick1Sp = thisClickSet(thisClickSet(:,2)==uSp(iSp),:); 
            [~,MTTidx,~] = intersect(MTT,thisClick1Sp(:,1));
            if length(unique(thisClick1Sp(:,2)))>1
                1;  % bad news, we didn't separate the species right
            end
            % can't load unevenly spaced events, so if they are uneven,
            % load one at a time for now. Could do something smarter, could
            % preallocate.
            if length(unique(diff(MTTidx)))>1
                thisSpecTemp = zeros(length(MTTidx),wMSP);
                thisWaveEnvTemp = zeros(length(MTTidx),wMSN);
                batches = split_into_batches(MTTidx);
                for iMTT = 1:length(batches)
                    thisBatch = batches{iMTT}(1):batches{iMTT}(2):batches{iMTT}(3);
                    [~,MTTidxloc,~] = intersect(MTTidx,thisBatch);
                    thisSpecTemp(MTTidxloc,:) = TPWSfileObj.MSP(thisBatch,:);
                    thisWaveEnvTemp(MTTidxloc,:) = TPWSfileObj.MSN(thisBatch,:);
                end
                thisSpec = mean(thisSpecTemp,1);
                thisWaveEnv = mean(abs(hilbert(thisWaveEnvTemp')),2)';
            else

                thisSpec = mean(TPWSfileObj.MSP(MTTidx,:),1);
                thisWaveEnv = mean(abs(hilbert(TPWSfileObj.MSN(MTTidx,:)')),2)';
            end
            
            [thisICI,~] = histc(diff(MTT(MTTidx))*24*60*60,ICIvec);
            if size(thisICI,2)>1
                thisICI = thisICI'; % wrong shape, happens when vector is all zeros.
            end
            allNewSpec = [allNewSpec;thisSpec];
            allNewWaveEnv = [allNewWaveEnv;thisWaveEnv];

            allNewICI = [allNewICI;thisICI'];
            allNewID = [allNewID;uSp(nSp)];
        end
            
               
    end
end