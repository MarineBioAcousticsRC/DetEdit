function detEdit(userFunc)

% detEdit.m

% Main script to display interface, it takes the user parameter settings
% and plots data in different panels to annotate the data.

% Input
%
%   userFunc - Script user parameter settings.
%       Optional, user wil be prompt to select a scrip if not provided 
%
% Examples:
%
% detEdit(@yourDataSettings)
% 
% detEdit

sizePoints = 9; % the current points
sizeBack = 5; % the background points
sizeFPR = 7; % the current points
colorPoints = [1 .84 0];

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

p = getParams(userFunc,'analysis','detEdit');

%% Define subfolder that fit specified iteration
if p.iterationNum > 1
    for id = 2: str2num(p.iterationNum) % iterate id times according to p.iterationNum
        subfolder = ['TPWS',num2str(id)];
        p.tpwsDir = (fullfile(p.tpwsDir,subfolder));
    end
end

%% Check if TPWS file exists
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
    error(sprintf('No files matching file prefix ''%s'' found!',detfn))
elseif length(fileMatchIdx)>1
    % if more than one match, throw error
    error(sprintf('Multiple TPWS files match the file prefix ''%s''.\n Make the prefix more specific.',detfn))
end

matchingFile = fileList{fileMatchIdx};

%% Handle Transfer Function
% add in transfer function if desired
if p.tfSelect > 0
    [tf,tffreq,tfuppc] = getTransfunc(p.filePrefix, p.tfName,p);
else
    tf = 0;
    disp('No TF Applied');
end
%% Generate FD, TD, and ID files if needed
[zFD,zID,fNameList] = buildLabelFiles(matchingFile, p.tpwsDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load detections and false detections
% MTT = time MPP = peak-peak % MSN = waveform %MSP = spectra
fNameList.TPWS = fullfile(p.tpwsDir,matchingFile);
load(fNameList.TPWS,'MTT','MPP')

% if you have more than "maxDetLoad" detections, don't load all spectra and
% time series into memory. You can access them from the disk instead.
% Note: Can't remove duplicates in that case, because matlab won't let you
% select non-contiguous sets from files stored on disk.
nDets = length(MTT);

if ~exist('maxDetLoad')
    loadMSP = true;
elseif isempty(maxDetLoad)
    loadMSP = true;
else
    loadMSP = nDets <= maxDetLoad; % true/false, if you have more detections than
end
% the maximum load this becomes false.
ic1 = [];
if loadMSP
    % remove duplicates from MTT (can't do this if too many detections to load into memory).
    [uMTT,ia1,ic1] = unique(MTT);
    if (length(uMTT) ~= length(MTT))
        disp([' TimeLevel Data NOT UNIQUE - removed:   ', ...
            num2str(length(ic1) - length(ia1))]);
    end
    load(fNameList.TPWS,'MSP','MSN')
else
    ia1 = [1:length(MTT)]';
end

[r,c] = size(MTT); %get shape of array
if (r > c)
    clickTimes = MTT(ia1);
    clickLevels = MPP(ia1);
else
    clickTimes = MTT(ia1)';
    clickLevels = MPP(ia1)';
end

if p.specploton && loadMSP
    % if p.specploton and there aren't too many detections, load spectra
    csn = MSN(ia1,:);
    csp = MSP(ia1,:);
else
    disp('Number of detections exceeds max for loading');
end

%% Apply tf and remove low amplitude detections
clickLevels = clickLevels + tf;
ib1 = find(clickLevels >= p.threshRL);

if (size(ib1,1) ~= size(clickLevels,1)) && ~loadMSP % catch for case where enforcing
    % min RL threshold on large dataset creates non-continuous indices.
    error('detEdit:RL',['Error: Re-run makeTPWS to enforce your minimum peak to peak RL threshold.\n',...
        'You cannot do it here because you have too many detections to load into memory.\n',...
        sprintf('TPWS minimum RL = %d \ndetEdit minimum RL = %d',min(clickLevels),p.threshRL)])
end

if (size(ib1,1) == 0)
    % min RL threshold excludes all detections
    error('detEdit:RL',['Error: No detections meet the minimum peak to peak RL threshold.\n',...
        sprintf('TPWS maximum RL = %d \ndetEdit minimum RL = %d',max(clickLevels),p.threshRL)])
end
% prune by RL only if spectra & waveforms have been loaded

if p.specploton && loadMSP
    disp([' Removed too low:',num2str(length(ia1)-length(ib1))]);
    clickTimes = clickTimes(ib1);
    clickLevels = clickLevels(ib1);
    keepers = ia1(ib1);
    
    csn = csn(ib1,:);
    csp = csp(ib1,:);
else
    keepers = ia1;
end

%% Make FD file intersect with MTT
load(fNameList.FD)  % false detection times zFD
zFD = rmUnmatchedDetections(MTT, zFD);
save(fNameList.FD,'zFD');

% Make ID file intersect with MTT
load(fNameList.ID)  % identified detection times zID
zID = rmUnmatchedDetections(MTT, zID);
save(fNameList.ID,'zID');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate bout starts and ends
% TODO: this is calculated in mkLTSA, we should save it there instead of
% recalculating here but would create backward compatibility issues

[nb,eb,sb,bd] = calculate_bouts(clickTimes,p);
bFlag = 0;

%% Make LTSA session file
lsfn = strrep(matchingFile,'TPWS','LTSA');
fNameList.LTSA = fullfile(p.tpwsDir,lsfn);
Altsa = exist(fNameList.LTSA,'file');
if Altsa ~= 2
    disp(['Error: LTSA Sessions File Does Not Exist: ',fNameList.LTSA])
    return
else
    disp('Loading LTSA Sessions, please wait ...')
    load(fNameList.LTSA)   % LTSA sessions: pwr and pt structures
    disp('Done Loading LTSA Sessions')
    sltsa = size(pt);
    if (sltsa(2) ~= nb)
        disp(['Error: Number of LTSA sessions calculated here doesn''t match ',...
            'input LTSA file. Check ltsaMax parameter.'])
        return
    end
end

%% Set up Tests for False Detections
% The false positive estimate tool picks every Nth click to test. If you
% have false positives in zFD, you can pick out only the true ones to
% determine which click indices to look at (that happens just below),
% but if the user then adds or removes anything from zFD, the indices won't
% adjust. Provide a warning to tell the user there might be an issue.
if ~isempty(zFD)
    disp(strcat('WARNING: This dataset contains false-flagged detections.  ', ...
        'Remove them using modDet prior to estimating false positive rate.'))
end

[~,trueClickIDx] = setdiff(clickTimes, zFD);
ixfd = (1: p.c4fd : length(trueClickIDx));  % selected to test for False Det
testClickIdx = trueClickIDx(ixfd);

A6 = exist(fNameList.TD,'file');
if (A6 ~= 2)
    zTD = -1.*ones(nb,4);
    save(fNameList.TD,'zTD');    % create new TD
    disp(' Make new TD file');
else
    load(fNameList.TD)
    if (length(zTD(:,1)) ~= nb)
        disp([' Problem with TD file:',fNameList.TD]);
        return
    end
end


%% Compute Spectra Plot Parameters
% max and min for LTSA frequency
fiminLTSA = 0;% TODO: make this configurable
fimaxLTSA = p.sampleRate/2 ; % in kHz 100 or 160 kHz

% set ltsa step size
iPwr = 1;
while isempty(pwr{1,iPwr}) && iPwr<length(pwr)
    iPwr = iPwr+1;
end

% ToDo: Seems like LTSA parameters (step size and frequency bins) could
% be calculated from LTSA directly, especially because there are
% inconsistent step sizes in some LTSAs.
if any(strcmp('dfManual',fieldnames(p)))&& ~isempty(p.dfManual)
    % allow non-standard ltsa step sizes
    df = p.dfManual;
else
    df = 1000*fimaxLTSA/(size(pwr{1,iPwr},1)-1);
end

% for LTSA PLOT
f = 1000*fiminLTSA:df:1000*fimaxLTSA;
if p.tfSelect > 0 % tfParams isn't being used...
    tfLTSA = interp1(tffreq,tfuppc,f,'linear','extrap')'; % add to LTSA vector
else
    tfLTSA = zeros(size(f))';
end

if p.specploton %Does anyone ever turn p.specploton off??
    % check length of MSP
    inFileMat = matfile(fNameList.TPWS);
    if ~loadMSP
        MSP = inFileMat.MSP(1,:);
    end
    smsp2 = size(MSP,2);% 2nd element is num fft points
    ift = 1:smsp2;
    % make frequency vector that matches spectral bins
    if any(strcmp('f',fieldnames(inFileMat)))
        fmsp = inFileMat.f;
    else
        fmsp = [];
    end
    if isempty(fmsp)
        fmsp = ((p.sampleRate/2)/(smsp2-1))*ift - (p.sampleRate/2)/(smsp2-1);
        fprintf('No freq vector in TPWS file. Using approximation based on sample rate.\n')
    end
    % find the indices that are in the range of interest
    fi = find(fmsp > p.fLow &...
        fmsp <= p.fHi);
    fimint = fi(1); fimaxt = fi(end);
    ft = fmsp(fi);
    
    % for the PP vs RMS plot
    if (p.tfSelect > 0)
        Ptfpp = interp1(tffreq,tfuppc,fmsp*1000,'linear','extrap');
    else
        Ptfpp = zeros(1,smsp2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = input('Starting Session:  ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause on
disp('Press ''b'' key to go backward ')
disp('Press any other key to go forward')
cc = ' ';  % avoids crash when first bout too short
yell = [];
blag = 0;
if (size(zTD,2) == 2) % this seems to patch on extra columns
    % to old zTD matrices that maybe only had the first two. Probably
    % only needed for backward compatibility
    zTD = [zTD,-1.*ones(length(zTD),2)];
    save(fNameList.TD,'zTD');
end

% Display ID legend:
% (note since it's not in the loop, if people close it,
% it won't come back in this session.
ID_legend(p)

%% Main Loop
% loop over the number of bouts (sessions)
onerun = 1; % What does this do?

while (k <= nb)
    disp([' BEGIN SESSION: ',num2str(k)]);
    % load in FD, ID and TD each session in case these have been modified
    load(fNameList.FD); % brings in zFD
    load(fNameList.ID); % brings in zID
    load(fNameList.TD); % brings in zTD
    
    % If all time series are loaded:
    % Make PP versus RMS plot for all clicks, if all time series are loaded
    figure(51); clf; set(51,'name',sprintf('RL pp vs. RL rms (left shift by %d)',p.threshRL))
    h51 = gca;
    
    % Make RMS versus frequency plot for all clicks
    figure(53); clf; set(53,'name',sprintf('RL rms vs. Peak freq. (left shift by %d)',p.threshRL))
    h53 = gca;
    
    if p.threshHiFreq ~=0 %any(strcmp('threshHiFreq',fieldnames(p)))
        ymax = p.threshHiFreq + 1;
    else
        ymax = fmsp(end);  % yaxis max of plot 53 (Default)
    end
    
    if p.specploton && loadMSP
        xmsp0All = csp + repmat(Ptfpp,size(csp,1),1);
        [xmspAll,im] = max(xmsp0All(:,fimint:fimaxt),[],2); % maximum between flow-100kHz
        
        % calculate peak-to-peak amplitude including transfer function
        if isrow(clickLevels)
            xmppAll = clickLevels-tf+ Ptfpp(im + fimint-1); % vectorized version
        else
            xmppAll = clickLevels'-tf+ Ptfpp(im + fimint-1); % vectorized version
        end
        % turn diagonal to vertical (easier way to find thresholds)
        if isrow(xmspAll)
            pxmspAll = xmspAll' - p.slope*(xmppAll - p.threshRL); %use slope of 1 to mod xmsp for plot
        elseif isrow(xmppAll)
            pxmspAll = xmspAll - p.slope*(xmppAll' - p.threshRL); %use slope of 1 to mod xmsp for plot
        else
            pxmspAll = xmspAll - p.slope*(xmppAll - p.threshRL); %use slope of 1 to mod xmsp for plot
        end
        plot(h51,pxmspAll,xmppAll,'o','MarkerSize',sizeBack,'MarkerEdgeColor',[.7,.7,.7],'UserData',clickTimes)
        title(h51,['Based on ',num2str(length(xmppAll)),' clicks']);
        
        % apply RMS threshold to figure (51)
        if (p.threshRMS > 0)
            if onerun == 1
                if p.threshPP > 0
                    badClickTime = clickTimes(pxmspAll < p.threshRMS &...
                        xmppAll' < p.threshPP);  % for all false if below RMS threshold
                else
                    badClickTime = clickTimes(pxmspAll < p.threshRMS);
                end
                disp(['Number of Detections Below RMS threshold = ',num2str(length(badClickTime))])
                zFD = [zFD; badClickTime];   % cummulative False Detection matrix
                save(fNameList.FD,'zFD')
            end
            if p.threshPP > 0
                xtline = [p.threshRMS,p.threshRMS]; ytline = [ min(xmppAll),p.threshPP];
            else
                xtline = [p.threshRMS,p.threshRMS]; ytline = [ min(xmppAll),max(xmppAll)];
            end
            hold(h51,'on');
            plot(h51,xtline,ytline,'r')
            hold(h51,'off');
            %p.threshRMS = 0;
        end
        
        % plot RMS vs frequency plot, keeping RMS vertical like in fig(51)
        freqAll = fmsp(im + fimint-1);
        plot(h53,pxmspAll,freqAll,'o','MarkerSize',sizeBack,'MarkerEdgeColor',[.7,.7,.7],'UserData',clickTimes)
        title(h53,['Based on total of ',num2str(length(freqAll)),' clicks']);
        % apply High Frequency threshold to figure (53)
        if onerun == 1
            if (p.threshHiFreq > 0)
                badClickTime = clickTimes(freqAll > p.threshHiFreq);  % for all false if below RMS threshold
                disp(['Number of Detections Below Freq threshold = ',num2str(length(badClickTime))])
                zFD = [zFD; badClickTime];   % cummulative False Detection matrix
                save(fNameList.FD,'zFD')
                %p.threshHiFreq = 0;
            end
        end
        if (p.threshHiFreq > 0)
            xtline = [min(pxmspAll),max(pxmspAll)]; ytline = [p.threshHiFreq ,p.threshHiFreq];
            hold(h53,'on');
            plot(h53,xtline,ytline,'r')
            hold(h53,'off');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find detections and false detections within this bout (session)
    J = []; JFD =[]; Jtrue = []; XFD = []; JID = [];
    J = find(clickTimes >= sb(k) & clickTimes <= eb(k));
    if p.specploton && loadMSP
        % have to load consecutive detections if reading from disk
        J = J(1):J(end);
        csnJ = csn(J,:);
        cspJ = csp(J,:);
    elseif p.specploton % only load the ones you need for this session
        csnJ = inFileMat.MSN(keepers(J),:);
        cspJ = inFileMat.MSP(keepers(J),:);
    end
    
    % get indices  of test clicks in this session
    XFD = find(clickTimes(testClickIdx) >= sb(k) &...
        clickTimes(testClickIdx) <= eb(k));
    zTD(k,1) = length(XFD);
    
    % Test for XFD and strcmp('x or z or w') - if no test points skip
    % x = true, z = false, w = window
    if (isempty(XFD) && (strcmp(cc,'x') || ...
            strcmp(cc,'z') || strcmp(cc,'w')))
        disp(' NO Test Detections, so skip')
        k = k + 1;
        continue
    end
    
    if ~isempty(J) % if there are detection in this session
        t = clickTimes(J); % detection times in this session
        disp([' Detection times:',num2str(length(t))]);
        if (~isempty(XFD))
            xt = clickTimes(testClickIdx(XFD));  %times to test for False Detection
            xPP = clickLevels(testClickIdx(XFD));   %amplitude for test False Detection
            disp([' Test False Detection times:',num2str(zTD(k,1))]),
        else
            xt = [];
        end
        RL = clickLevels(J);         % received levels in this session
        nd = length(J);     % number of detections in this session
        
        % get false detection times that intersect with detection times
        K2 = []; % holds false indices
        ff2 = 0;
        tfd = [];
        if (~isempty(zFD)) % get times and indices of false detections
            [tfd,K2,~] = intersect(t,zFD(:,1));
            rlFD = RL(K2);
        end
        if ~isempty(K2) % if this session contains false detections
            ff2 = 1; % set false flag to true
            if p.specploton
                wavFD = norm_wav(mean(csnJ(K2,:),1)); % calculate mean false time series
                specFD = cspJ(K2,:); % get set of false spectra
            end
            disp([' False detections:',num2str(length(K2))])
        else
            ff2 = 0;
            disp(' No False Detections')
        end
        
        % get ID'd detection times that intersect with detection times
        K3 = []; % holds Id'd indices
        ff3 = 0; % becomes positive if you have ID's detections in this session
        tID = []; % times of ID'd detections
        IDidx = [];
        if ~isempty(zID)
            [tID,K3,IDidx] = intersect(t,zID(:,1));
            rlID = RL(K3);
        end
        if ~isempty(K3)
            ff3 = 1;
            spCodeSet = zID(IDidx,2); % get ID codes for everything in this session
            specIDs = unique(spCodeSet); % get unique ID codes
            if p.specploton % get mean spectra for each ID'd type
                wavID = [];
                specID = [];
                for iSpID = 1:length(specIDs)
                    thisSet = spCodeSet == specIDs(iSpID);
                    wavID(iSpID,:) = norm_wav(mean(csnJ(K3(thisSet,:),:),1));
                    specID(iSpID,:) = mean(cspJ(K3(thisSet,:),:),1);
                end
            end
            disp([' ID detections:',num2str(length(K3))])
        else
            ff3 = 0;
            disp(' No identified detections (ID)')
        end
        
        
        % Calculate indices of detections which are neither false, ID'd,
        JFD = J(K2);
        JID = J(K3);
        JFIM = union(JFD,JID);
        [Jtrue,iJ,~]= setxor(J,JFIM); % find all true detections
        trueTimes = clickTimes(Jtrue);% vector of true times in this session
        
        if p.specploton
            cspJtrue = cspJ(iJ,:); % true spectra in this session
            csnJtrue = csnJ(iJ,:); % true time series in this session
            wtrue = norm_wav(nanmean(csnJtrue,1)); % mean of true spectra in this session
            strue = nanmean(cspJtrue,1); % mean of true time series in this session
        end
        disp([' True Detections: ',num2str(length(trueTimes))])
    else
        disp('Error: no detections between bout start and end')
        return
    end
    dt = diff(t)*24*60*60; % inter-detection interval (IDI) and convert from days to seconds
    
    if ff2 % calculate IDI for false and id'd detections
        dtFD = dt(K2(1:end-1));
    end
    if ff3
        dtID = dt(K3(1:end-1));
    end
    disp(['END SESSION: ',num2str(k),' Start: ',datestr(sb(k)),...
        ' End:',datestr(eb(k))])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate number of detections per bin
    [KB,binCX,binT,binC] = ndets_per_bin(t,xt,RL,dt,p.minNdet,nd,p.binDur);
    % filter empty and low number bins
    if isempty(KB) % not sure what this case does?
        disp(['No bins with at least ',num2str(p.minNdet),' detections'])
        binT = 0;
        binRL = 0;
        binC = 0;
        k = k + 1;  % go to next
        continue
    end
    if (strcmp(cc,'w') && (zTD(k,2) == 0))
        disp(['Session: ',num2str(k),' # Test Detect Bins: ',...
            num2str(length(binCX)),' but NO False']);
        zTD(k,3) = length(binCX);
        zTD(k,4) = 0;
        save(fNameList.TD,'zTD');
        k = k + 1;
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number detection per spectral bin in LTSA
    % make a spectra in figure 50
    PT = pt{1,k};   % LTSA session time vector
    pwr1 = pwr{1,k};  % LTSA power vector
    nbinS = length(PT);
    if (nbinS == 0)
        disp('No LTSA for this Session');
        PT(1) = sb(k) ; PT(2) = eb(k); % start end times for plots
        pwr1(1:length(f)) = ones; % make uniform LTSA
    else
        pwr1 = pwr1((1000*fiminLTSA/df)+1:round(1000*fimaxLTSA/df)+1,:);
    end
    durS = PT(end) - PT(1);
    
    if p.specploton
        % allSPEC = norm_spec(cspJ,fimint,fimint,fimaxt);
        figure(50);clf;set(50,'name','Frequency Spectra')
        h50 = gca;
        figure(52);clf;set(52,'name','Waveform')
        h52 = gca;
        trueSpec = [];
        if ~isempty(trueTimes)
            % plot average true click spectrum
            trueSpec = norm_spec_simple(cspJtrue,fimint,fimaxt);
            plot(h50,ft,trueSpec,'Linewidth',4)
            % average true click waveform
            plot(h52, wtrue);
        else
            disp(['No true with at least ',num2str(p.minNdet),' detections'])
        end
        if ff2   % average false click spec
            SPEC2 = norm_spec_simple(specFD,fimint,fimaxt);
            % plot average false click spectrum
            hold(h50, 'on')
            plot(h50,ft,SPEC2,'r','Linewidth',4)
            hold(h50, 'off')
            % plot average false click waveform
            hold(h52, 'on')
            plot(h52,wavFD + 0.5 ,'r');
            hold(h52, 'off')
        end
        if ff3  % average id click spec
            specID_norm = [];
            for iSpec = 1:size(specID,1)
                specID_norm(iSpec,:) = norm_spec_simple(specID(iSpec,:),fimint,fimaxt);
            end
            % plot average ID'd click spectra
            hold(h50, 'on')
            hID = plot(h50,ft,specID_norm,'Linewidth',4);
            hold(h50, 'off')
            
            % plot average ID'd click waveform(s)
            hold(h52, 'on')
            hID2 = plot(h52,(wavID + repmat(-1*rand(size(hID)),1,length(wavID)))');
            
            for iC = 1:length(hID) % set colors
                set(hID(iC),'Color',p.colorTab(specIDs(iC),:))
                set(hID2(iC),'Color',p.colorTab(specIDs(iC),:))
            end
            hold(h52, 'off')
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add detections of this session in figure 51 and 53
        % for all detections in this session, calculate xmpp and xmsp
        xmsp0 = cspJ + repmat(Ptfpp,size(cspJ,1),1); % add transfer fun to session's spectra
        [xmsp,im] = max(xmsp0(:,fimint:fimaxt),[],2);
        if isrow(RL)
            xmpp = RL - tf + Ptfpp([im + fimint - 1]);
        else
            xmpp = RL' - tf + Ptfpp([im + fimint - 1]);
        end
        % turn diagonal to vertical
        if ~isempty(xmsp) && ~isempty(xmpp)
            if isrow(xmsp)
                pxmsp = xmsp' - p.slope*(xmpp - p.threshRL); %use slope of 1 to mod xmsp for plot
            elseif isrow(xmpp)
                pxmsp = xmsp - p.slope*(xmpp' - p.threshRL);
            else
                pxmsp = xmsp - p.slope*(xmpp - p.threshRL);
            end
        end
        
        % Plot  PP versus RMS Plot for this session
        hold(h51, 'on')
        h51.ColorOrderIndex = 1;
        plot(h51,pxmsp,xmpp,'.','MarkerSize',sizePoints,'UserData',t)% true ones in blue
        if ~loadMSP % plot threshold line now because no background data
            if (p.threshRMS > 0)
                if onerun == 1
                    if p.threshPP > 0
                        badClickTime = t(pxmsp < p.threshRMS &...
                            xmpp' < p.threshPP);  % for all false if below RMS threshold
                    else
                        badClickTime = t(pxmsp < p.threshRMS);
                    end
                    disp(['Number of Detections Below RMS threshold = ',num2str(length(badClickTime))])
                    zFD = [zFD; badClickTime];   % cummulative False Detection matrix
                    save(fNameList.FD,'zFD')
                    if ~isempty(zFD) % get times and indices of false detections
                        [tfd,K2,~] = intersect(t,zFD(:,1));
                        rlFD = RL(K2);
                    end
                    if ~isempty(K2) % if this session contains false detections
                        ff2 = 1; % set false flag to true
                        if p.specploton
                            wavFD = norm_wav(mean(csnJ(K2,:),1)); % calculate mean false time series
                            specFD = cspJ(K2,:); % get set of false spectra
                        end
                        dtFD = dt(K2(1:end-1));
                        disp([' False detections:',num2str(length(K2))])
                    else
                        ff2 = 0;
                        disp(' No False Detections')
                    end
                end
            end
            if p.threshPP > 0 && exist('plotaxes','var')
                xtline = [p.threshRMS,p.threshRMS]; ytline = [ plotaxes.minPP,p.threshPP];
            elseif p.threshPP > 0
                xtline = [p.threshRMS,p.threshRMS]; ytline = [ min(xmpp),p.threshPP];
            else
                xtline = [p.threshRMS,p.threshRMS]; ytline = [ min(xmpp),max(xmpp)];
            end
            plot(h51,xtline,ytline,'r')
        end
        
        if ff2 % false in red
            plot(h51,pxmsp(K2),xmpp(K2),'r.','MarkerSize',sizePoints,'UserData',t(K2))
        end
        if ff3 % ID'd in associated color
            for iC2 = 1:length(specIDs) % set colors
                thisIDset = spCodeSet ==specIDs(iC2);
                hPP = plot(h51,pxmsp(K3(thisIDset)),xmpp(K3(thisIDset)),'.','MarkerSize',sizePoints,'UserData',t(K3(thisIDset)));
                set(hPP,'Color',p.colorTab(specIDs(iC2),:))
            end
        end
        
        hold(h51, 'off')
        
        % Plot RMS vs frequency plot for this session
        hold(h53, 'on')
        h53.ColorOrderIndex = 1;
        freq = fmsp(im + fimint -1);
        plot(h53,pxmsp,freq,'.','MarkerSize',sizePoints,'UserData',t) % true ones in blue
        if ~loadMSP
            if onerun == 1
                if (p.threshHiFreq > 0)
                    badClickTime = t(freq > p.threshHiFreq);  % for all false if below RMS threshold
                    disp(['Number of Detections Below Freq threshold = ',num2str(length(badClickTime))])
                    zFD = [zFD; badClickTime];   % cummulative False Detection matrix
                    save(fNameList.FD,'zFD')
                    if ~isempty(zFD) % get times and indices of false detections
                        [tfd,K2,~] = intersect(t,zFD(:,1));
                        rlFD = RL(K2);
                    end
                    if ~isempty(K2) % if this session contains false detections
                        ff2 = 1; % set false flag to true
                        if p.specploton
                            wavFD = norm_wav(mean(csnJ(K2,:),1)); % calculate mean false time series
                            specFD = cspJ(K2,:); % get set of false spectra
                        end
                        disp([' False detections:',num2str(length(K2))])
                        dtFD = dt(K2(1:end-1));
                    else
                        ff2 = 0;
                        disp(' No False Detections')
                    end
                end
            end
            if p.threshHiFreq > 0 && exist('plotaxes','var')
                xtline = [plotaxes.minRMS,plotaxes.maxRMS]; ytline = [p.threshHiFreq ,p.threshHiFreq];
            elseif p.threshHiFreq > 0
                xtline = [min(pxmsp),max(pxmsp)]; ytline = [p.threshHiFreq ,p.threshHiFreq];
            end
            plot(h53,xtline,ytline,'r')
        end
        if ff2 % false in red
            plot(h53,pxmsp(K2),freq(K2),'r.','MarkerSize',sizePoints,'UserData',t(K2))
        end
        if ff3 % ID'd in associated color
            for iC2 = 1:length(specIDs) % set colors
                thisIDset = spCodeSet ==specIDs(iC2);
                hPP = plot(h53,pxmsp(K3(thisIDset)),freq(K3(thisIDset)),'.','MarkerSize',sizePoints,'UserData',t(K3(thisIDset)));
                set(hPP,'Color',p.colorTab(specIDs(iC2),:))
            end
        end

        hold(h53, 'off')
        if p.threshHiFreq > 0
            ylim(h53,[p.fLow ymax+10])
        else
            ylim(h53,[p.fLow ymax])
        end
    end
    
    % add figure labels
    xlabel(h50,'Frequency (kHz)');
    grid(h50,'on')
    xlim(h50, 'manual');
    ylim(h50,[0 1]);
    xlim(h50,[p.fLow,p.fHi])
    
    xlabel(h51,'dB RMS')
    ylabel(h51,'dB Peak-to-peak')
    if exist('plotaxes','var')
        xlim(h51,[plotaxes.minRMS,plotaxes.maxRMS])
        ylim(h51,[plotaxes.minPP,plotaxes.maxPP])
    end
    
    xlabel(h52,'Time (1ms @ 200kHz)');
    ylabel(h52,' Normalized Amplitude');
    
    xlabel(h53,'dB RMS')
    ylabel(h53,'Peak Frequency (kHz)')
    if exist('plotaxes','var')
        xlim(h53,[plotaxes.minRMS,plotaxes.maxRMS])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots stuff now in figure(201)
    warning('off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(201);clf
    colormap(201, jet)
    
    set(201,'name','LTSA and time series')
    hA201 = subplot_layout; % Top panel, Figure 201: Received Level
    plot(hA201(1),t,RL,'.','MarkerSize',sizePoints,'UserData',t)
    hold(hA201(1),'on')
    if ff2 % plot False detections in red
        plot(hA201(1),tfd,rlFD,'r.','MarkerSize',sizePoints,'UserData',tfd)
        % disp([' false det plotted:',num2str(length(tfd))])
    end
    if ff3 % plot ID'd detections in associated color
        spCodeSet = zID(IDidx,2); % get species codes for everything in this session
        specIDs = unique(spCodeSet); % get unigue species codes
        for iC2 = 1:length(specIDs) % set colors
            thisIDset = spCodeSet ==specIDs(iC2);
            hRLID = plot(hA201(1),tID(thisIDset),rlID(thisIDset),'.',...
                'MarkerSize',sizePoints,'UserData',tID(thisIDset));
            set(hRLID,'Color',p.colorTab(specIDs(iC2),:))
        end
    end
    hold(hA201(1),'off')
    axis(hA201(1),[PT(1) PT(end) p.rlLow p.rlHi])
    datetick(hA201(1),'x',15,'keeplimits')
    grid(hA201(1),'on')
    tstr(1) = {fNameList.TPWS};
    tstr(2) = {['Session: ',num2str(k),'/',num2str(nb),' Start Time ',...
        datestr(sb(k)),' Detect = ',num2str(nd)]};
    title(hA201(1),tstr);
    ylabel(hA201(1),'RL [dB re 1\muPa]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % middle panel LTSA
    c = (p.ltsaContrast/100) .* pwr1 + p.ltsaBright;
    image(PT,f/1000,c,'parent',hA201(2))
    set(hA201(2),'yDir','normal')
    axis(hA201(2),[PT(1) PT(end) p.ltsaLims])%v2(4)
    ylabel(hA201(2),'Frequency (kHz)')
    datetick(hA201(2),'keeplimits')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bottom panel, Figure 201: Inter-Detection Interval
    % make two copies of dt points for brush
    tdt2 = [];
    dt2 = [];
    ldt = length(dt);
    if ldt > 0
        tdt2 = reshape([t(1:ldt),t((1:ldt)+1)]',2*ldt,1);
        dt2 = reshape([dt,dt]',2*ldt,1);
        
        hICI = plot(hA201(3),tdt2,dt2,'.');
        set(hICI,'MarkerSize',sizePoints,'MarkerFaceColor','b','LineStyle','none','UserData',tdt2)
        % Note: plotyy is buggy in 2012b, axis handles work only if called
        % using "axes" and avoid calls to "subplot"
        
        % Do setup for 1st axes
        axis(hA201(3),[PT(1) PT(end) 0 p.dtHi])
        datetick(hA201(3),'x',15,'keeplimits')
        Ytick = 0:p.dtHi/10:p.dtHi;
        set(hA201(3),'YTick',Ytick)
        datetick(hA201(3),'x',15,'keeplimits')
        grid(hA201(3),'on')
        ylabel(hA201(3),'Time between detections [s]')
%         title(AX(1),'Inter-Detection Interval (IDI)')
        
        %%% plot FD, ID
        hold(hA201(3),'on')
        if ff2
            plot(hA201(3),tfd(2:end),dtFD,'.r','MarkerSize',sizePoints,'UserData',tfd(2:end))
            % no need to double FD since only the blue points are brush captured
        end
        if ff3 % plot ID'd in associated color
            for iC3 = 1:length(specIDs) % set colors
                thisIDset = spCodeSet ==specIDs(iC3);
                hICI = plot(hA201(3),tID(thisIDset(2:end)),...
                    dtID(thisIDset(2:end)),'.','MarkerSize',sizePoints,...
                    'UserData',tID(thisIDset));
                set(hICI,'Color',p.colorTab(specIDs(iC2),:))
            end
        end
        hold(hA201(3),'off')
    else
        plot(0,0);
    end
    
    % if you have items brushed in yellow, highlight those on each plot
    if p.specploton && ~isempty(yell) && ~isempty(csnJ)
        hold(hA201(1),'on')
        plot(hA201(1),t(yell),RL(yell),'ko','MarkerSize',sizeBack,'UserData',t(yell));
        hold(hA201(1),'off');
        
        % for diffs, yell can't exceed length dt, which could happen if you
        % grabbed the last point in the vector, so:
        yellDT = yell(yell<length(dt));
        hold(hA201(3),'on')
        plot(hA201(3),t(yellDT),dt(yellDT),'ko','MarkerSize',sizeBack,'UserData',t(yell));
        hold(hA201(3),'off')
        
        hold(h50,'on') % add click to spec plot in BLACK
        cspJy = mean(cspJ(yell,:),1);
        tSPEC = norm_spec_simple(cspJy,fimint,fimaxt);
        plot(h50,ft,tSPEC,'k','Linewidth',4);
        hold(h50,'off')
        
        hold(h51,'on')
        plot(h51,pxmsp(yell),xmpp(yell),'ko','MarkerSize',sizeBack,...
            'LineWidth',2,'UserData',clickTimes(K2))
        hold(h51,'off')
        
        hold(h52,'on') % add click to waveform plot in BLACK
        plot(h52,norm_wav(mean(csnJ(yell,:),1))' + 1.5,'k');
        hold(h52,'off')
        
        hold(h53,'on')
        plot(h53,pxmsp(yell),freq(yell),'ko','MarkerSize',sizeBack,...
            'LineWidth',2,'UserData',clickTimes(K2))
        hold(h53,'off')
    end
    % end of plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yell = [];
    pause  % wait for user input.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    onerun = onerun+1;
    % get key stroke
    cc = get(gcf,'CurrentCharacter');
    
    % if brush selected get key
    if strcmp(cc,'p')
       h = brush;
       set(h,'Color',[1,1,0],'Enable','on'); % light yellow [.9290 .6940 .1250]
       waitfor(gcf,'CurrentCharacter')
       set(h,'Enable','off')
       cc = get(gcf,'CurrentCharacter');
    end
    
    if strcmp(cc,'u') || strcmp(cc,'g') || strcmp(cc,'r')
        % detections were flagged by user
        disp(' Update Display') % Stay on same bout
        % get brushed data and figure out what to do based on color:
        [yell,zFD,zID,bFlag] = brush_color(gca,cc,zFD,zID,p.colorTab,t);
        
    elseif strcmp(cc,'s') % change time diff scale on bottom plot of 201
        p.dtHi = input(' Update IPI scale (sec):  '); % Set IPI scale
        
    elseif strcmp(cc,'d') % change RL scale on top plot of 201
        p.rlLow = input(' Update RL low (dB):  '); % Set RL low
        p.rlHi = input(' Update RL high (dB):  '); % Set RL high
        
    elseif strcmp(cc,'a')% change LTSA parameters
        p.ltsaContrast = input(sprintf('  Current Contrast %d. Update Contrast:  ',p.ltsaContrast));
        p.ltsaBright = input(sprintf('  Current Brightness %d. Update Brightness:  ',p.ltsaBright));
        
    elseif strcmp(cc,'<') % change RMS threshold plot 51
        p.threshRMS = input(' Set RMS Threshold:  '); % Set false thres
        
    elseif strcmp(cc,':') % change RMS threshold plot 51
        p.threshPP = input(' Set PP Threshold:  '); % Set false thres
        
    elseif strcmp(cc,'^') % change High Frequency threshold plot 51
        p.threshHiFreq = input(' Set High Frequency Threshold:  '); % Set false thres
        
    elseif strcmp(cc,'!') % change High Frequency threshold plot 51
        ymax = input(' Update High Frequency axis:  '); % Set false thres
        
    elseif strcmp(cc,'b') % Move backward one bout
        if k ~= 1
            k = k-1;
        end
        onerun = 1;
    elseif strcmp(cc,'f') % assign ALL as false
        disp(['Number of False Detections Added = ',num2str(length(trueTimes))])
        if ~isempty(zID)
            %[newFD,~] = setdiff(t,zID(:,1)); % remove from zID
            [~,iCID] = setdiff(zID(:,1),t); % remove from zID
            zID = zID(iCID,:);
        end
        newFD = t;
        zFD = [zFD; newFD; trueTimes]; % Add everything to zFD
        
    elseif strcmp(cc,'t') %assign ALL as true
        [zFD,~] = setdiff(zFD(:,1),t);
        disp(['Remaining False Detections = ',num2str(length(zFD))])
        if ~isempty(zID)
            [~,iCID] = setdiff(zID(:,1),t);
            zID = zID(iCID,:);
        end
    elseif strcmp(cc,'j')% jump to non-consecutive session
        prompt = 'Jump to Session: ';
        kjump = input(prompt);
        if (kjump > 0 && kjump < nb)
            k = kjump;
        end
        onerun = 1;
        
    elseif (strcmp(cc,'x') || strcmp(cc,'z') ) % test click for random False Detect
        if ~isempty(XFD)
            zTD(k,2) = 0;
            for inxfd = 1 : zTD(k,1)
                hold(hA201(1),'on')
                testTimes = xt(inxfd);
                xH201a = plot(hA201(1),testTimes,xPP(inxfd),'o','MarkerEdgeColor',colorPoints,...
                    'MarkerSize',sizeFPR,'LineWidth',2);
                hold(hA201(1),'off')
                inxfdDT = inxfd(inxfd<length(dt));
                hold(hA201(3),'on')
                xH201b = plot(hA201(3),testTimes,dt(inxfdDT),'o','MarkerEdgeColor',colorPoints,...
                    'MarkerSize',sizeFPR,'LineWidth',2);
                hold(hA201(3),'off')
                disp(['Showing #: ',num2str(inxfd),' click. Press ''z'' to reject']);
                if (p.specploton == 1)
                    hold(h50,'on')  % add click to spec plot in BLACK
                    %plot(h50,ft,trueSpec,'Linewidth',2);
                    clickInBoutIdx = find(t==testTimes);
                    testSnip = csnJtrue(clickInBoutIdx,:);
                    testSpectrum = cspJtrue(clickInBoutIdx,:);
                    
                    
                    % make low freq part = 0
                    tempSPEC = norm_spec_simple(testSpectrum,fimint,fimaxt);
                    xH0 = plot(h50,ft,tempSPEC,'Color',colorPoints,'Linewidth',4);
                    hold(h50,'off')
                    
                    hold(h52,'on') % add click to waveform plot in BLACK
                    xH2 = plot(h52,norm_wav(testSnip)' + 1.5,'Color',colorPoints,...
                        'Linewidth',2);
                    hold(h52,'off')
                    
                    hold(h51,'on')
                    % get click index relative to bout
                    xH1 = plot(h51,pxmsp(clickInBoutIdx),xmpp(clickInBoutIdx),...
                        'o','MarkerEdgeColor',colorPoints,'MarkerSize',sizeFPR,...
                        'LineWidth',2);
                    hold(h51,'off')
                    
                    hold(h53,'on')
                    xH3 = plot(h53,pxmsp(clickInBoutIdx),freq(clickInBoutIdx),...
                        'o','MarkerEdgeColor',colorPoints,'MarkerSize',sizeFPR,...
                        'LineWidth',2);
                    hold(h53,'off')
                end
                pause
                cc = get(gcf,'CurrentCharacter');
                % fill evaluated points
                xH201a.MarkerFaceColor = colorPoints;
                xH201b.MarkerFaceColor = colorPoints;
                xH1.MarkerFaceColor = colorPoints;
                xH3.MarkerFaceColor = colorPoints;
                if (strcmp(cc,'z'))
                    zTD(k,2) = zTD(k,2) + 1;
                    zFD = [zFD; xt(inxfd)]; % add to FD
                end
                delete([xH0,xH1,xH2,xH3])
            end
            disp([' Tested: ',num2str(zTD(k,1)),' False: ',...
                num2str(zTD(k,2))]);
            
        end
        k = k+1;
    elseif (strcmp(cc,'w') && (zTD(k,2) > 0))  % test 5 min window
        % Test 5 min window
        zTD = test_false_bins(k,zTD,xt,xPP,binCX);
        k = k+1;
        
    else
        k = k+1;  % move forward one bout
        onerun = 1;
    end
    
    % after edits, remove duplicate labels and save updated vectors
    if ~isempty(zFD)
        zFD = unique(zFD);
    end
    if ~isempty(zID)
        [~,uniqueID] = unique(zID(:,1));
        zID = zID(uniqueID,:);
    end
    save(fNameList.FD,'zFD')
    save(fNameList.ID,'zID')
    save(fNameList.TD,'zTD')
    
    % don't end if you used paintbrush on last record
    if (k > nb) && bFlag
        k = nb;
        disp(' Last Record')
        
    end
    bFlag = 0;
end
pause off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make zFD unique
uZFD = [];  ia = []; ic = [];
load(fNameList.FD);   % load false detections
[uzFD,ia,ic] = unique(zFD);     % make zFD have unique entries
if (length(ia) ~= length(ic))
    disp([' False Detect NOT UNIQUE - removed:   ', ...
        num2str(length(ic) - length(ia))]);
end
zFD = uzFD;
save(fNameList.FD,'zFD');
tfinal = find(zTD(:,1) > 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(['Number of Starting Detections = ',num2str(length(clickTimes)+2)])
disp(' ')
disp(['Number of True Detections = ',num2str(length(clickTimes)-length(zFD)+2)])
disp(' ')
disp(['Number of False Detections = ',num2str(length(zFD)-1)])
disp(' ')
disp(['Number of Test Detections & False Detect = ',num2str(sum(zTD(tfinal,:)))])
disp(' ')
disp(['Done with file ',fNameList.TPWS])

commandwindow;