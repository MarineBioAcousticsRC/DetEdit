function compute_transformations

global p dPARAMS

if p.loadMSP % if all spectra are loaded into memory
    % compute peak freq directly
    [~,im] = max(dPARAMS.csp(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2); % maximum between flow-100kHz
    
    % compute rms
    cspLinear = 10.^(dPARAMS.csp/10);
    binWidth = (dPARAMS.ft(2)-dPARAMS.ft(1));%Fs/nfft;
    dPARAMS.RMSall = 10*log10(sum(cspLinear(:,dPARAMS.fimint:dPARAMS.fimaxt)...
        .*binWidth,2)) + dPARAMS.tf;
    
else
    % otherwise, load spectra in batches and calculate peak freq
    batchList = 0:p.maxDetLoad:length(dPARAMS.clickTimes);
    if batchList(end) ~= length(dPARAMS.clickTimes)
        % make sure last number makes it into list
        batchList = [batchList,length(dPARAMS.clickTimes)];
    end
    dPARAMS.xmspAll = [];
    dPARAMS.RMSall = [];
    im = [];
    for iBatch = 1:length(batchList)-1
        thisBatch = (batchList(iBatch)+1):batchList(iBatch+1);
        cspTemp = dPARAMS.inFileMat.MSP(thisBatch,:);
        [~,imTemp] = max(cspTemp(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2);
        
        % compute rms
        cspTempLinear = 10.^(cspTemp/10);
        binWidth = (dPARAMS.ft(2)-dPARAMS.ft(1));%Fs/nfft;
        RMSallTemp = 10*log10(sum(cspTempLinear(:,dPARAMS.fimint:dPARAMS.fimaxt)...
            .*binWidth,2)) + dPARAMS.tf;
        
        dPARAMS.RMSall = [dPARAMS.RMSall;RMSallTemp];
        im = [im;imTemp];
    end
end
dPARAMS.xmppAll = dPARAMS.clickLevels - dPARAMS.tf+ ...
    dPARAMS.Ptfpp(im + dPARAMS.fimint-1)'; % vectorized version
dPARAMS.transfRMSall = dPARAMS.RMSall - p.slope*(dPARAMS.xmppAll - p.threshRL); %use slope of 1 to mod xmsp for plot

dPARAMS.freqAll = dPARAMS.fmsp(im + dPARAMS.fimint-1);
