function boutMotion

global dPARAMS p fNameList zID zTD zFD dHANDLES


disp([' BEGIN SESSION: ',num2str(dPARAMS.k)]);
% load in FD, ID and TD each session in case these have been modified
load(fNameList.FD); % brings in zFD
load(fNameList.ID); % brings in zID
load(fNameList.TD); % brings in zTD



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find detections and false detections within this session
J = []; JFD =[]; Jtrue = []; dPARAMS.XFD = []; JID = [];
J = find(dPARAMS.clickTimes >= dPARAMS.sb(dPARAMS.k) & dPARAMS.clickTimes <= dPARAMS.eb(dPARAMS.k));

if p.loadMSP
    % have to load consecutive detections if reading from disk
    J = J(1):J(end);
    dPARAMS.csnJ = dPARAMS.csn(J,:);
    dPARAMS.cspJ = dPARAMS.csp(J,:);
else % only load the ones you need for this session
    dPARAMS.csnJ = dPARAMS.inFileMat.MSN(dPARAMS.keepers(J),:);
    dPARAMS.cspJ = dPARAMS.inFileMat.MSP(dPARAMS.keepers(J),:);
end

% get indices  of test clicks in this session
dPARAMS.XFD = find(dPARAMS.clickTimes(dPARAMS.testClickIdx) >= dPARAMS.sb(dPARAMS.k) &...
    dPARAMS.clickTimes(dPARAMS.testClickIdx) <= dPARAMS.eb(dPARAMS.k));
zTD(dPARAMS.k,1) = length(dPARAMS.XFD);
save(fNameList.TD,'zTD')

% Test for XFD and strcmp('x or z or w') - if no test points skip
% x = true, z = false, w = window
if (isempty(dPARAMS.XFD) && (strcmp(dPARAMS.cc,'x') || ...
        strcmp(dPARAMS.cc,'z') || strcmp(dPARAMS.cc,'w')))
    disp(' NO Test Detections, so skip')
    dPARAMS.k = dPARAMS.k + 1;
    return
end

if ~isempty(J) % if there are detection in this session
    dPARAMS.t = dPARAMS.clickTimes(J); % detection times in this session
    disp([' Detection times:',num2str(length(dPARAMS.t))]);
    if (~isempty(dPARAMS.XFD))
        dPARAMS.xt = dPARAMS.clickTimes(dPARAMS.testClickIdx(dPARAMS.XFD));  %times to test for False Detection
        dPARAMS.xPP = dPARAMS.clickLevels(dPARAMS.testClickIdx(dPARAMS.XFD));   %amplitude for test False Detection
        disp([' Test False Detection times:',num2str(zTD(dPARAMS.k,1))]),
    else
        dPARAMS.xt = [];
    end
    dPARAMS.RL = dPARAMS.clickLevels(J);         % received levels in this session
    dPARAMS.nd = length(J);     % number of detections in this session
    
    % get false detection times that intersect with detection times
    dPARAMS.K2 = []; % holds false indices
    dPARAMS.ff2 = 0;
    dPARAMS.tfd = [];
    if (~isempty(zFD)) % get times and indices of false detections
        [dPARAMS.tfd,dPARAMS.K2,~] = intersect(dPARAMS.t,zFD(:,1));
        dPARAMS.rlFD = dPARAMS.RL(dPARAMS.K2);
    end
    if ~isempty(dPARAMS.K2) % if this session contains false detections
        dPARAMS.ff2 = 1; % set false flag to true
        if p.specploton
            dPARAMS.wavFD =  norm_wav(mean(dPARAMS.csnJ(dPARAMS.K2,:),1)); % calculate mean false time series
            dPARAMS.specFD = dPARAMS.cspJ(dPARAMS.K2,:); % get set of false spectra
        end
        disp([' False detections:',num2str(length(dPARAMS.K2))])
    else
        dPARAMS.ff2 = 0;
        disp(' No False Detections')
    end
    
    % get ID'd detection times that intersect with detection times
    dPARAMS.K3 = []; % holds Id'd indices
    dPARAMS.ff3 = 0; % becomes positive if you have ID's detections in this session
    dPARAMS.tID = []; % times of ID'd detections
    dPARAMS.labelConf = []; % label confidence scores
    IDidx = [];
    if ~isempty(zID)
        [dPARAMS.tID,dPARAMS.K3,IDidx] = intersect(dPARAMS.t,zID(:,1));
        dPARAMS.rlID = dPARAMS.RL(dPARAMS.K3);
        if size(zID,2)>2
            % store label confidence if available, otherwise fill with na.
            dPARAMS.labelConf = round(zID(IDidx,3)*1000)/1000;
        else
            dPARAMS.labelConf = nan(length(IDidx),1);
        end
        dPARAMS.labelConfIdx = (dPARAMS.labelConf>=p.minLabelConfidence |...
            isnan(dPARAMS.labelConf));
        dPARAMS.unlabeledIdx = 1:length(dPARAMS.t);
        dPARAMS.unlabeledIdx(dPARAMS.K3(dPARAMS.labelConfIdx)) = [];
    else
        dPARAMS.unlabeledIdx = 1:length(dPARAMS.t);
    end
    if ~isempty(dPARAMS.K3)
        dPARAMS.ff3 = 1;
        dPARAMS.spCodeSet = zID(IDidx,2); % get ID codes for everything in this session
        dPARAMS.specIDs = unique(dPARAMS.spCodeSet); % get unique ID codes

        disp([' Labels above confidence limit:',num2str(sum(dPARAMS.labelConfIdx))])

        if p.specploton % get mean spectra for each ID'd type
            dPARAMS.wavID = [];
            dPARAMS.specID = [];
            goodIdx = find(dPARAMS.labelConfIdx);
            for iSpID = 1:length(dPARAMS.specIDs)
                dPARAMS.thisIDset{iSpID} = intersect(find(dPARAMS.spCodeSet == dPARAMS.specIDs(iSpID)),goodIdx);
                dPARAMS.wavID(iSpID,:) =  norm_wav(mean(dPARAMS.csnJ(dPARAMS.K3(dPARAMS.thisIDset{iSpID},:),:),1));
                dPARAMS.specID(iSpID,:) = mean(dPARAMS.cspJ(dPARAMS.K3(dPARAMS.thisIDset{iSpID} ,:),:),1);
            end
        end
        disp([' ID detections:',num2str(length(dPARAMS.K3))])
    else
        dPARAMS.ff3 = 0;
        disp(' No identified detections (ID)')
    end
    
    % Calculate indices of detections which are neither false nor ID'd
    JFD = J(dPARAMS.K2);
    if isfield(dPARAMS,'labelConfIdx')
        JID = J(dPARAMS.K3(dPARAMS.labelConfIdx));
    else
        JID = J(dPARAMS.K3);
    end
    JFIM = union(JFD,JID);
    [Jtrue,iJ,~]= setxor(J,JFIM); % find all true detections
    dPARAMS.trueTimes = dPARAMS.clickTimes(Jtrue);% vector of true times in this session
    
    dPARAMS.cspJtrue = dPARAMS.cspJ(iJ,:); % true spectra in this session
    dPARAMS.csnJtrue = dPARAMS.csnJ(iJ,:); % true time series in this session
    dPARAMS.wtrue = norm_wav(nanmean(dPARAMS.csnJtrue,1)); % mean of true spectra in this session
    dPARAMS.strue = nanmean(dPARAMS.cspJtrue,1); % mean of true time series in this session
    
    disp([' True Detections: ',num2str(length(dPARAMS.trueTimes))])
else
    disp('Error: no detections between bout start and end')
    return
end
dPARAMS.dt = diff(dPARAMS.t)*24*60*60; % inter-detection interval (IDI) and convert from days to seconds
dPARAMS.dtUnlabeled = diff(dPARAMS.t(dPARAMS.unlabeledIdx))*24*60*60;

if dPARAMS.ff2 % calculate IDI for false and id'd detections
    dPARAMS.dtFD = dPARAMS.dt(dPARAMS.K2(1:end-1));
end
if dPARAMS.ff3
    dPARAMS.dtID = dPARAMS.dt(dPARAMS.K3(1:end-1));
end

disp(['END SESSION: ',num2str(dPARAMS.k),' Start: ',datestr(dPARAMS.sb(dPARAMS.k)),...
    ' End:',datestr(dPARAMS.eb(dPARAMS.k))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate number of detections per bin
[KB,binCX,binT,binC] = ndets_per_bin(dPARAMS.t,dPARAMS.xt,dPARAMS.RL,dPARAMS.dt,p.minNdet,dPARAMS.nd,p.binDur);
% filter empty and low number bins
if isempty(KB) % not sure what this case does?
    disp(['No bins with at least ',num2str(p.minNdet),' detections'])
    binT = 0;
    binRL = 0;
    binC = 0;
    dPARAMS.k = dPARAMS.k + 1;  % go to next
    return
end
if (strcmp(dPARAMS.cc,'w') && (zTD(dPARAMS.k,2) == 0))
    disp(['Session: ',num2str(dPARAMS.k),' # Test Detect Bins: ',...
        num2str(length(binCX)),' but NO False']);
    zTD(dPARAMS.k,3) = length(binCX);
    zTD(dPARAMS.k,4) = 0;
    save(fNameList.TD,'zTD');
    dPARAMS.k = dPARAMS.k + 1;
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number detection per spectral bin in LTSA
% make a spectra in figure 50
dPARAMS.PT(1) = dPARAMS.sb(dPARAMS.k) ;
dPARAMS.PT(2) = dPARAMS.eb(dPARAMS.k); % start end times for plots
dPARAMS.pwr1 = dPARAMS.pwr{1,dPARAMS.k};  % LTSA power vector
nbinS = length(dPARAMS.pwr1 );
if (nbinS == 0)
    disp('No LTSA for this Session');
    dPARAMS.pwr1(1:length(dPARAMS.f)) = ones; % make uniform LTSA
% else
%     dPARAMS.pwr1 = dPARAMS.pwr1((1000*dPARAMS.fiminLTSA/dPARAMS.df)+1:round(1000*dPARAMS.fimaxLTSA/dPARAMS.df)+1,:);
end
dPARAMS.durS = dPARAMS.PT(end) - dPARAMS.PT(1);


dPARAMS.trueSpec = [];
if ~isempty(dPARAMS.trueTimes)
    % plot average true click spectrum
    dPARAMS.trueSpec = norm_spec_simple(dPARAMS.cspJtrue,dPARAMS.fimint,dPARAMS.fimaxt);
else
    disp(['No true with at least ',num2str(p.minNdet),' detections'])
end
if dPARAMS.ff2   % average false click spec
    dPARAMS.SPEC2 = norm_spec_simple(dPARAMS.specFD,dPARAMS.fimint,dPARAMS.fimaxt);
end
if dPARAMS.ff3  % average id click spec
    % dPARAMS.spCodeSet = zID(IDidx,2); % get species codes for everything in this session
    dPARAMS.specIDs = unique(dPARAMS.spCodeSet); % get unigue species codes
    dPARAMS.specID_norm = [];
    
    for iSpec = 1:size(dPARAMS.specIDs,1)
        dPARAMS.specID_norm(iSpec,:) = norm_spec_simple(dPARAMS.specID(iSpec,:),...
            dPARAMS.fimint,dPARAMS.fimaxt);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add detections of this session in figure 51 and 53
% for all detections in this session, calculate xmpp and xmsp
%%% Does this have to be calculated every time? %%%
xmsp0 = dPARAMS.cspJ + repmat(dPARAMS.Ptfpp,size(dPARAMS.cspJ,1),1); % add transfer fun to session's spectra
[~,im] = max(dPARAMS.cspJ(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2);
cspJLinear = 10.^(dPARAMS.cspJ/10);
binWidth = (dPARAMS.ft(2)-dPARAMS.ft(1));%Fs/nfft;
RMS = 10*log10(sum(cspJLinear(:,dPARAMS.fimint:dPARAMS.fimaxt).*binWidth,2))+ dPARAMS.tf; % maximum between flow-100kHz
dPARAMS.xmpp = dPARAMS.RL' - dPARAMS.tf + dPARAMS.Ptfpp([im + dPARAMS.fimint - 1]);
[xmsp,im] = max(xmsp0(:,dPARAMS.fimint:dPARAMS.fimaxt),[],2);
dPARAMS.pxmsp = xmsp - p.slope*(dPARAMS.xmpp' - p.threshRL);


dPARAMS.transfRMS = RMS - p.slope*(dPARAMS.xmpp' - p.threshRL);
%%%-----%%%

if ~p.loadMSP % plot threshold line now because no background data
    if (p.threshRMS > 0)
        
        if ~isempty(zFD) % get times and indices of false detections
            [dPARAMS.tfd,dPARAMS.K2,~] = intersect(dPARAMS.t,zFD(:,1));
            dPARAMS.rlFD = dPARAMS.RL(dPARAMS.K2);
        end
        if ~isempty(dPARAMS.K2) % if this session contains false detections
            dPARAMS.ff2 = 1; %%% these should be set only once, not in figure
            if p.specploton
                dPARAMS.wavFD = norm_wav(mean(dPARAMS.csnJ(dPARAMS.K2,:),1)); % calculate mean false time series
                dPARAMS.specFD = dPARAMS.cspJ(dPARAMS.K2,:); % get set of false spectra
            end
            dPARAMS.dtFD = dPARAMS.dt(dPARAMS.K2(1:end-1));
            disp([' False detections:',num2str(length(dPARAMS.K2))])
        else
            dPARAMS.ff2 = 0;%%% these should be set only once, not in figure
            disp(' No False Detections')
        end
    end
end


dPARAMS.freq = dPARAMS.fmsp(im + dPARAMS.fimint -1);
if ~p.loadMSP
    if (p.threshHiFreq > 0)
        if ~isempty(zFD) % get times and indices of false detections
            [dPARAMS.tfd,dPARAMS.K2,~] = intersect(dPARAMS.t,zFD(:,1));
            dPARAMS.rlFD = dPARAMS.RL(dPARAMS.K2);
        end
        if ~isempty(dPARAMS.K2) % if this session contains false detections
            dPARAMS.ff2 = 1; % set false flag to true
            dPARAMS.wavFD =  norm_wav(mean(dPARAMS.csnJ(dPARAMS.K2,:),1)); % calculate mean false time series
            dPARAMS.specFD = dPARAMS.cspJ(dPARAMS.K2,:); % get set of false spectra
            
            disp([' False detections:',num2str(length(dPARAMS.K2))])
            dPARAMS.dtFD = dPARAMS.dt(dPARAMS.K2(1:end-1));
        else
            dPARAMS.ff2 = 0;
            disp(' No False Detections')
        end
    end
end

figure201

figure51
figure53

figure50
figure52

if isfield(dHANDLES,'hID')
    % bring legend to top
    figure(dHANDLES.hID)
end
