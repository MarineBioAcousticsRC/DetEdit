function apply_auto_thresh

global p dPARAMS zFD fNameList
% verify autothresh

applyCode = input('Confirm applying frequency and RMS thresholds to entire file (y/n):',...
    's');

if strcmpi(applyCode,'y')
    disp('Applying automatic thresholds')
    
    % apply RMS threshold
    % calculate peak-to-peak amplitude including transfer function
    if (p.threshHiFreq > 0)
        badClickTime = dPARAMS.clickTimes(dPARAMS.freqAll > p.threshHiFreq);  % for all false if below freq threshold
        disp(['Number of Detections Below Freq threshold = ',num2str(length(badClickTime))])
        zFD = [zFD; badClickTime];   % cummulative False Detection matrix
    end
    
    if p.threshPP > 0
        badClickTime = dPARAMS.clickTimes(dPARAMS.pxmspAll < p.threshRMS &...
            dPARAMS.xmppAll < p.threshPP);  % for all false if below RMS threshold
    else
        badClickTime = dPARAMS.clickTimes(dPARAMS.pxmspAll < p.threshRMS);
    end
    disp(['Number of Detections Below RMS threshold = ',num2str(length(badClickTime))])
    zFD = [zFD; badClickTime];   % cummulative False Detection matrix
    zFD = unique(zFD);
    save(fNameList.FD,'zFD')
else
    disp('No thresholds applied')
end
