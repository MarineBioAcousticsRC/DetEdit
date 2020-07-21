function figure53


global dPARAMS dHANDLES p

% Make RMS versus frequency plot for all clicks
figure(dHANDLES.RMSvFreqfig);clf
dHANDLES.h53 = gca;


% plot RMS vs frequency plot, keeping RMS vertical like in fig(51)
if p.loadMSP
    plot(dHANDLES.h53,dPARAMS.transfRMSall,dPARAMS.freqAll,'o',...
        'MarkerSize',p.sizeBack,'MarkerEdgeColor',[.7,.7,.7],...
        'UserData',dPARAMS.clickTimes)
    title(dHANDLES.h53,['Based on total of ',num2str(length(dPARAMS.freqAll)),' clicks']);
    % apply High Frequency threshold to figure (53)
    if (p.threshHiFreq > 0)
        xtline = [min(dPARAMS.transfRMSall),max(dPARAMS.transfRMSall)];
        ytline = [p.threshHiFreq ,p.threshHiFreq];
        hold(dHANDLES.h53,'on');
        plot(dHANDLES.h53,xtline,ytline,'r')
        hold(dHANDLES.h53,'off');
    end
end
% Plot RMS vs frequency plot for this session
dHANDLES.h53.ColorOrderIndex = 1;
hold(dHANDLES.h53, 'on')
dHANDLES.RMS53 = plot(dHANDLES.h53,dPARAMS.transfRMS(dPARAMS.unlabeledIdx),...
    dPARAMS.freq(dPARAMS.unlabeledIdx),'.','MarkerSize',...
    p.sizePoints,'UserData',dPARAMS.t(dPARAMS.unlabeledIdx),'Visible',dPARAMS.NoLabel_Toggle); % true ones in blue

if ~p.loadMSP
    if p.threshHiFreq > 0 %&& isfield(dPARAMS,'plotaxes')
%         xtline = [p.rmsLow,p.rmsHi];
%         ytline = [p.threshHiFreq ,p.threshHiFreq];
%     elseif p.threshHiFreq > 0
        xtline = [min(dPARAMS.transfRMS),max(dPARAMS.transfRMS)];
        ytline = [p.threshHiFreq ,p.threshHiFreq];
    
        plot(dHANDLES.h53,xtline,ytline,'r')
    end
end
dHANDLES.RMSFD53 = [];
if dPARAMS.ff2 % false in red
    dHANDLES.RMSFD53 = plot(dHANDLES.h53,dPARAMS.transfRMS(dPARAMS.K2),dPARAMS.freq(dPARAMS.K2),...
        'r.','MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K2),...
        'Visible',dPARAMS.FD_Toggle);
end

dHANDLES.RMSID53 =  cell(size(p.colorTab,1),1);
if dPARAMS.ff3 % ID'd in associated color
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        iColor = dPARAMS.specIDs(iC2);
        thisIDset = dPARAMS.spCodeSet == dPARAMS.specIDs(iC2);
        dHANDLES.RMSID53{iColor} = plot(dHANDLES.h53,dPARAMS.transfRMS(dPARAMS.K3(thisIDset)),...
            dPARAMS.freq(dPARAMS.K3(thisIDset)),'.',...
            'MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K3(thisIDset)));
        set(dHANDLES.RMSID53{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC2),:),...
            'Visible',dPARAMS.ID_Toggle{iColor});
    end
end


hold(dHANDLES.h53, 'off')
if p.threshHiFreq > 0
    ylim(dHANDLES.h53,[p.fLow p.threshHiFreq+10])
else
    ylim(dHANDLES.h53,[p.fLow p.fHi])
end
xlim(dHANDLES.h53,[p.rmsLow,p.rmsHi])

xlabel(dHANDLES.h53,'Transformed Received Level (dB_r_m_s re 1\muPa)')
ylabel(dHANDLES.h53,'Peak Frequency (kHz)')

if ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h53,'on')
    plot(dHANDLES.h53,dPARAMS.transfRMS(dPARAMS.yell),dPARAMS.freq(dPARAMS.yell),'ko',...
        'MarkerSize',p.sizeBack,...
        'LineWidth',2,'UserData',dPARAMS.clickTimes(dPARAMS.K2))
    hold(dHANDLES.h53,'off')
end