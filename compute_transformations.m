function compute_transformations

global p dPARAMS

if p.loadMSP % if all spectra are loaded into memory
    % compute peak freq directly
    dPARAMS.xmsp0All = dPARAMS.csp + repmat(dPARAMS.Ptfpp,size(dPARAMS.csp,1),1);
    %     [dPARAMS.xmspAll,im] = max(dPARAMS.xmsp0All(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2); % maximum between flow-100kHz
    [~,im] = max(dPARAMS.xmsp0All(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2); % maximum between flow-100kHz
    
    linearsp = 10.^(dPARAMS.xmsp0All/20);
    binWidth = (dPARAMS.ft(2)-dPARAMS.ft(1))*1000;%Fs/nfft;
    dPARAMS.xmspAll = 10*log10(sum(linearsp(:,dPARAMS.fimint:dPARAMS.fimaxt).*binWidth,2)); % maximum between flow-100kHz
    
else
    % otherwise, load spectra in batches and calculate peak freq
    batchList = 0:p.maxDetLoad:length(dPARAMS.clickTimes);
    if batchList(end) ~= length(dPARAMS.clickTimes)
        % make sure last number makes it into list
        batchList = [batchList,length(dPARAMS.clickTimes)];
    end
    dPARAMS.xmspAll = [];
    im = [];
    for iBatch = 1:length(batchList)-1
        thisBatch = (batchList(iBatch)+1):batchList(iBatch+1);
        cspTemp = dPARAMS.inFileMat.MSP(thisBatch,:);
        xmsp0AllTemp = cspTemp + repmat(dPARAMS.Ptfpp,size(cspTemp,1),1);
        %         [xmspAllTemp,imTemp] = max(xmsp0AllTemp(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2);
        [~,imTemp] = max(xmsp0AllTemp(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2);
        
        linearsp = 10.^(dPARAMS.xmsp0AllTemp/20);
        binWidth = (dPARAMS.ft(2)-dPARAMS.ft(1))*1000;%Fs/nfft;
        dPARAMS.xmspAllTemp = 10*log10(sum(linearsp(:,dPARAMS.fimint:dPARAMS.fimaxt).*binWidth,2)); % maximum between flow-100kHz
        
        dPARAMS.xmspAll = [dPARAMS.xmspAll;xmspAllTemp];
        im = [im;imTemp];
    end
end
dPARAMS.xmppAll = dPARAMS.clickLevels - dPARAMS.tf+ ...
    dPARAMS.Ptfpp(im + dPARAMS.fimint-1)'; % vectorized version
dPARAMS.pxmspAll = dPARAMS.xmspAll - p.slope*(dPARAMS.clickLevels - p.threshRL); %use slope of 1 to mod xmsp for plot

dPARAMS.freqAll = dPARAMS.fmsp(im + dPARAMS.fimint-1);
