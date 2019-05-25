function keyAction(varargin)

global zID zFD zTD fNameList p dPARAMS

disp('Key action detected')
dPARAMS.cc = get(gcf,'CurrentCharacter');
% if brush selected get key
% if strcmp(dPARAMS.cc,'p')
%     h = brush;
%     set(h,'Color',[1,1,0],'Enable','on'); % light yellow [.9290 .6940 .1250]
%     waitfor(gcf,'CurrentCharacter')
%     set(h,'Enable','off')
%     dPARAMS.cc = get(gcf,'CurrentCharacter');
% end

if strcmp(dPARAMS.cc,'u') || strcmp(dPARAMS.cc,'g') || strcmp(dPARAMS.cc,'r')
    % detections were flagged by user
    disp(' Update Display') % Stay on same bout
    % get brushed data and figure out what to do based on color:
    [dPARAMS.yell,zFD,zID,dPARAMS.bFlag] = brush_color(gca,dPARAMS.cc,zFD,zID,p.colorTab,dPARAMS.t);
    
elseif ~isempty(str2num(dPARAMS.cc))
    
    [zFD,zID] = brush_color_number(gca,dPARAMS.cc,zFD,zID,p.colorTab,dPARAMS.t);
    
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
%    if ~isempty(dPARAMS.XFD)
%         zTD(k,2) = 0;
%         for inxfd = 1 : zTD(k,1)
%             hold(hA201(1),'on')
%             testTimes = xt(inxfd);
%             xH201a = plot(hA201(1),testTimes,xPP(inxfd),'o','MarkerEdgeColor',colorPoints,...
%                 'MarkerSize',sizeFPR,'LineWidth',2);
%             hold(hA201(1),'off')
%             inxfdDT = inxfd(inxfd<length(dt));
%             hold(hA201(3),'on')
%             xH201b = plot(hA201(3),testTimes,dt(inxfdDT),'o','MarkerEdgeColor',colorPoints,...
%                 'MarkerSize',sizeFPR,'LineWidth',2);
%             hold(hA201(3),'off')
%             disp(['Showing #: ',num2str(inxfd),' click. Press ''z'' to reject']);
%             if (p.specploton == 1)
%                 hold(h50,'on')  % add click to spec plot in BLACK
%                 %plot(h50,ft,trueSpec,'Linewidth',2);
%                 clickInBoutIdx = find(t==testTimes);
%                 testSnip = csnJtrue(clickInBoutIdx,:);
%                 testSpectrum = cspJtrue(clickInBoutIdx,:);
%                 
%                 
%                 % make low freq part = 0
%                 tempSPEC = norm_spec_simple(testSpectrum,fimint,fimaxt);
%                 xH0 = plot(h50,ft,tempSPEC,'Color',colorPoints,'Linewidth',4);
%                 hold(h50,'off')
%                 
%                 hold(h52,'on') % add click to waveform plot in BLACK
%                 xH2 = plot(h52,norm_wav(testSnip)' + 1.5,'Color',colorPoints,...
%                     'Linewidth',2);
%                 hold(h52,'off')
%                 
%                 hold(h51,'on')
%                 % get click index relative to bout
%                 xH1 = plot(h51,pxmsp(clickInBoutIdx),xmpp(clickInBoutIdx),...
%                     'o','MarkerEdgeColor',colorPoints,'MarkerSize',sizeFPR,...
%                     'LineWidth',2);
%                 hold(h51,'off')
%                 
%                 hold(h53,'on')
%                 xH3 = plot(h53,pxmsp(clickInBoutIdx),freq(clickInBoutIdx),...
%                     'o','MarkerEdgeColor',colorPoints,'MarkerSize',sizeFPR,...
%                     'LineWidth',2);
%                 hold(h53,'off')
%             end
%             pause
%             dPARAMS.cc = get(gcf,'CurrentCharacter');
%             % fill evaluated points
%             xH201a.MarkerFaceColor = colorPoints;
%             xH201b.MarkerFaceColor = colorPoints;
%             xH1.MarkerFaceColor = colorPoints;
%             xH3.MarkerFaceColor = colorPoints;
%             if (strcmp(dPARAMS.cc,'z'))
%                 zTD(k,2) = zTD(k,2) + 1;
%                 zFD = [zFD; xt(inxfd)]; % add to FD
%             end
%             delete([xH0,xH1,xH2,xH3])
%         end
%         disp([' Tested: ',num2str(zTD(k,1)),' False: ',...
%             num2str(zTD(k,2))]);
%         
%     end
    dPARAMS.k = dPARAMS.k+1;
elseif (strcmp(dPARAMS.cc,'w') && (zTD(k,2) > 0))  % test 5 min window
    % Test 5 min window
    zTD = test_false_bins(k,zTD,dPARAMS.xt,dPARAMS.xPP,dPARAMS.binCX);
    k = k+1;
    
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
        k = k+1;  % move forward one bout
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