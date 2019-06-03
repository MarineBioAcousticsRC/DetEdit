function keyAction(varargin)

global zID zFD zTD fNameList p dPARAMS dHANDLES

disp('Key action detected')
dPARAMS.cc = get(gcf,'CurrentCharacter');
dPARAMS.yell = [];
% if brush selected get key
% if strcmp(dPARAMS.cc,'p')
%     h = brush;
%     set(h,'Color',[1,1,0],'Enable','on'); % light yellow [.9290 .6940 .1250]
%     waitfor(gcf,'CurrentCharacter')
%     set(h,'Enable','off')
%     dPARAMS.cc = get(gcf,'CurrentCharacter');
% end

if strcmp(dPARAMS.cc,'u') || strcmp(dPARAMS.cc,'g') || strcmp(dPARAMS.cc,'r')...
        || strcmp(dPARAMS.cc,'y')
    % detections were flagged by user
    disp(' Update Display') % Stay on same bout
    % get brushed data and figure out what to do based on color:
    [dPARAMS.yell,zFD,zID,dPARAMS.bFlag] = brush_color(gca,dPARAMS.cc,zFD,zID,p.colorTab,dPARAMS.t);
    
elseif ~isempty(str2num(dPARAMS.cc))
    
    [zFD,zID] = brush_color_number(gca,str2num(dPARAMS.cc),zFD,zID,p.colorTab,dPARAMS.t);
    
elseif strcmp(dPARAMS.cc,'s') % change time diff scale on bottom plot of 201
    p.dtHi = input(' Update IPI scale (sec):  '); % Set IPI scale
    
elseif strcmp(dPARAMS.cc,'d') % change RL scale on top plot of 201
    p.rlLow = input(' Update RL low (dB):  '); % Set RL low
    p.rlHi = input(' Update RL high (dB):  '); % Set RL high
    
elseif strcmp(dPARAMS.cc,'a')% change LTSA parameters
    p.ltsaContrast = input(sprintf('  Current Contrast %d. Update Contrast:  ',p.ltsaContrast));
    p.ltsaBright = input(sprintf('  Current Brightness %d. Update Brightness:  ',p.ltsaBright));
    
elseif strcmp(dPARAMS.cc,'<') % change RMS threshold plot 51
    p.threshRMS = input(' Set RMS Threshold:  '); % Set false thres
    
elseif strcmp(dPARAMS.cc,':') % change RMS threshold plot 51
    p.threshPP = input(' Set PP Threshold:  '); % Set false thres
    
elseif strcmp(dPARAMS.cc,'^') % change High Frequency threshold plot 51
    p.threshHiFreq = input(' Set High Frequency Threshold:  '); % Set false thres
    
elseif strcmp(dPARAMS.cc,'!') % change High Frequency threshold plot 51
    dPARAMS.ymax = input(' Update High Frequency axis:  '); % Set false thres
    
elseif strcmp(dPARAMS.cc,'b') % Move backward one bout
    if dPARAMS.k ~= 1
        dPARAMS.k = dPARAMS.k-1;
    end
    dPARAMS.onerun = 1;
elseif strcmp(dPARAMS.cc,'f') % assign ALL as false
    disp(['Number of False Detections Added = ',num2str(length(dPARAMS.trueTimes))])
    if ~isempty(zID)
        %[newFD,~] = setdiff(t,zID(:,1)); % remove from zID
        [~,iCID] = setdiff(zID(:,1),dPARAMS.t); % remove from zID
        zID = zID(iCID,:);
    end
    newFD = dPARAMS.t;
    zFD = [zFD; newFD; dPARAMS.trueTimes]; % Add everything to zFD
    
elseif strcmp(dPARAMS.cc,'t') %assign ALL as true
    [zFD,~] = setdiff(zFD(:,1),dPARAMS.t);
    disp(['Remaining False Detections = ',num2str(length(zFD))])
    if ~isempty(zID)
        [~,iCID] = setdiff(zID(:,1),dPARAMS.t);
        zID = zID(iCID,:);
    end
elseif strcmp(dPARAMS.cc,'j')% jump to non-consecutive session
    prompt = 'Jump to Session: ';
    kjump = input(prompt);
    if (kjump > 0 && kjump < dPARAMS.nb)
        dPARAMS.k = kjump;
    end
    dPARAMS.onerun = 1;
    
elseif (strcmp(dPARAMS.cc,'x') || strcmp(dPARAMS.cc,'z') ) % test click for random False Detect
    if ~isempty(dPARAMS.XFD)
        zTD(dPARAMS.k,2) = 0;
        % entering test mode, temporarily disable other callbacks to avoid
        % confusion with keys enered here
        set(dHANDLES.LTSAfig, 'KeyPressFcn',[])
        set(dHANDLES.RMSvPPfig, 'KeyPressFcn',[])
        set(dHANDLES.RMSvFreqfig, 'KeyPressFcn',[])
        set(dHANDLES.spectrafig, 'KeyPressFcn',[])
        set(dHANDLES.wavefig, 'KeyPressFcn',[])
        for inxfd = 1:zTD(dPARAMS.k,1)
            % update LTSA to show highlighted click
            hold(dHANDLES.LTSAsubs(1),'on')
            testTimes = dPARAMS.xt(inxfd);
            xH201a = plot(dHANDLES.LTSAsubs(1),testTimes,dPARAMS.xPP(inxfd),'o','MarkerEdgeColor',...
                p.colorPoints,'MarkerSize',p.sizeFPR,'LineWidth',2);
            hold(dHANDLES.LTSAsubs(1),'off')
            inxfdDT = inxfd(inxfd<length(dPARAMS.dt));
            hold(dHANDLES.LTSAsubs(3),'on')
            xH201b = plot(dHANDLES.LTSAsubs(3),testTimes,dPARAMS.dt(inxfdDT),...
                'o','MarkerEdgeColor',p.colorPoints,'MarkerSize',p.sizeFPR,'LineWidth',2);
            hold(dHANDLES.LTSAsubs(3),'off')
            disp(['Showing #: ',num2str(inxfd),' click. Press ''z'' to reject']);
            
            % update spectra
            hold(dHANDLES.h50,'on')  % add click to spec plot in BLACK
            plot(dHANDLES.h50,dPARAMS.ft,dPARAMS.trueSpec,'Linewidth',2);
            clickInBoutIdx = find(dPARAMS.trueTimes==testTimes);
            testSnip = dPARAMS.csnJtrue(clickInBoutIdx,:);
            testSpectrum = dPARAMS.cspJtrue(clickInBoutIdx,:);
            
            % Normalize spectrum
            tempSPEC = norm_spec_simple(testSpectrum,dPARAMS.fimint,dPARAMS.fimaxt);
            xH0 = plot(dHANDLES.h50,dPARAMS.ft,tempSPEC,'Color',p.colorPoints,'Linewidth',4);
            hold(dHANDLES.h50,'off')
            
            % Update waveform
            hold(dHANDLES.h52,'on') % add click to waveform plot in BLACK
            xH2 = plot(dHANDLES.h52,norm_wav(testSnip)' + 1.5,'Color',p.colorPoints,...
                'Linewidth',2);
            hold(dHANDLES.h52,'off')
            
            % Update rms v pp
            hold(dHANDLES.h51,'on')
            % get click index relative to bout
            xH1 = plot(dHANDLES.h51,dPARAMS.pxmsp(clickInBoutIdx),dPARAMS.xmpp(clickInBoutIdx),...
                'o','MarkerEdgeColor',p.colorPoints,'MarkerSize',p.sizeFPR,...
                'LineWidth',2);
            hold(dHANDLES.h51,'off')
            
            % Update rms v freq.
            hold(dHANDLES.h53,'on')
            xH3 = plot(dHANDLES.h53,dPARAMS.pxmsp(clickInBoutIdx),dPARAMS.freq(clickInBoutIdx),...
                'o','MarkerEdgeColor',p.colorPoints,'MarkerSize',p.sizeFPR,...
                'LineWidth',2);
            hold(dHANDLES.h53,'off')
            
            % wait for key
            pause
            
            dPARAMS.cc = get(gcf,'CurrentCharacter');
            
            % fill evaluated points
            xH201a.MarkerFaceColor = p.colorPoints;
            xH201b.MarkerFaceColor = p.colorPoints;
            xH1.MarkerFaceColor = p.colorPoints;
            xH3.MarkerFaceColor = p.colorPoints;
            
            if (strcmp(dPARAMS.cc,'z'))
                zTD(k,2) = zTD(k,2) + 1;
                zFD = [zFD; xt(inxfd)]; % add to FD
            end
            
            delete([xH0,xH1,xH2,xH3]) % can we just overwrite and delete after loop?
        end
        disp([' Tested: ',num2str(zTD(dPARAMS.k,1)),' False: ',...
            num2str(zTD(dPARAMS.k,2))]);
        
    end
    dPARAMS.k = dPARAMS.k+1;
    % re-eneable normal callbacks
    set(dHANDLES.LTSAfig, 'KeyPressFcn',@keyAction)
    set(dHANDLES.RMSvPPfig, 'KeyPressFcn',@keyAction)
    set(dHANDLES.RMSvFreqfig, 'KeyPressFcn',@keyAction)
    set(dHANDLES.spectrafig, 'KeyPressFcn',@keyAction)
    set(dHANDLES.wavefig, 'KeyPressFcn',@keyAction)
elseif (strcmp(dPARAMS.cc,'w') && (zTD(k,2) > 0))  % test 5 min window
    % Test 5 min window
    zTD = test_false_bins(k,zTD,dPARAMS.xt,dPARAMS.xPP,dPARAMS.binCX);
    k = k+1;
    
elseif strcmp(dPARAMS.cc,'e') % re-code one species ID with another
    % detect if data have been brushed, otherwise use whole set.
    [brushDate, ~, ~] = get_brushed(gca);
    if isempty(brushDate)
        tEdit = dPARAMS.t;
    else
        tTemp = brushDate;
        tEdit = intersect(tTemp,dPARAMS.t);
    end
    oldID = input(' Enter the ID you want to overwrite:  ');
    newID  = input(' Enter the ID you want to change it to (enter 0 for no ID):  ');
    addFlag = 0; % flag gets turned to 1 if we have to append to zID rather than change existing IDs
    if oldID == 0 %get everything that's unlabeled
        addFlag = 1;
        [dates2Append,~] = setdiff(tEdit,[dPARAMS.tfd;dPARAMS.tID]);
        
    elseif oldID ==99 % get everything that's false
        addFlag = 1;
        [dates2Append,iCFD] = intersect(zFD(:,1),tEdit);
        zFD(iCFD) = [];
    else
        [~,iCID] = intersect(zID(:,1),tEdit);
        oldIDLocs = find(zID(iCID,2)==oldID);
    end
    
    if newID == 0  && ~addFlag % user wants previously ID'd thing changed to unlabeled
        zID(iCID(oldIDLocs),:) = [];
    elseif newID == 99 && ~addFlag % user wants previously ID'd thing changed to false
        zFD = [zFD;zID(iCID(oldIDLocs),1)];
        zID(iCID(oldIDLocs),:) = [];
    elseif newID == 99 && oldID == 0 % user wants to change unlabeled to false
        zFD = [zFD;dates2Append];
    elseif newID ==0 && oldID == 99 % user wants to change false to unlabeled
        % nothing needs to happen, because they've already been removed
        % from zFD above.
    else
        if addFlag % user wants previously false or unlabeled thing changed to ID
            zID = [zID; [dates2Append,ones(size(dates2Append))*newID]];
        else  % user wants previously ID'dthing changed to different ID
            zID(iCID(oldIDLocs),2) = newID;
        end
    end
    
    % check for accidental zeros in ID, this can happen if people enter a
    % 0, remove those from the ID set.
    accidentalZeros = zID(:,2)==0;
    zID(accidentalZeros,:)=[];
else
    if dPARAMS.k == dPARAMS.nb
        uZFD = [];  ia = []; ic = [];
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
        disp(['Number of Starting Detections = ',num2str(length(dPARAMS.clickTimes)+2)])
        disp(' ')
        disp(['Number of True Detections = ',num2str(length(dPARAMS.clickTimes)-length(zFD)+2)])
        disp(' ')
        disp(['Number of False Detections = ',num2str(length(zFD)-1)])
        disp(' ')
        disp(['Number of Test Detections & False Detect = ',num2str(sum(zTD(tfinal,:)))])
        disp(' ')
        disp(['Done with file ',fNameList.TPWS])
        commandwindow
        return
    else
        dPARAMS.k = dPARAMS.k+1;  % move forward one bout
        onerun = 1;
    end
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
if (dPARAMS.k > dPARAMS.nb) && dPARAMS.bFlag
    dPARAMS.k = dPARAMS.nb;
    disp(' Last Record')
    
end
dPARAMS.bFlag = 0;
boutMotion