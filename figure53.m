function figure53


global dPARAMS dHANDLES p

% Make RMS versus frequency plot for all clicks
figure(dHANDLES.RMSvFreqfig);clf
dHANDLES.h53 = gca;


% plot RMS vs frequency plot, keeping RMS vertical like in fig(51)
if p.loadMSP
    plot(dHANDLES.h53,dPARAMS.pxmspAll,dPARAMS.freqAll,'o',...
        'MarkerSize',p.sizeBack,'MarkerEdgeColor',[.7,.7,.7],'UserData',dPARAMS.clickTimes)
    title(dHANDLES.h53,['Based on total of ',num2str(length(dPARAMS.freqAll)),' clicks']);
    % apply High Frequency threshold to figure (53)
    if (p.threshHiFreq > 0)
        xtline = [min(dPARAMS.pxmspAll),max(dPARAMS.pxmspAll)];
        ytline = [p.threshHiFreq ,p.threshHiFreq];
        hold(dHANDLES.h53,'on');
        plot(dHANDLES.h53,xtline,ytline,'r')
        hold(dHANDLES.h53,'off');
    end
end
% Plot RMS vs frequency plot for this session
dHANDLES.h53.ColorOrderIndex = 1;
hold(dHANDLES.h53, 'on')
plot(dHANDLES.h53,dPARAMS.pxmsp,dPARAMS.freq,'.','MarkerSize',...
    p.sizePoints,'UserData',dPARAMS.t) % true ones in blue

if ~p.loadMSP
    if p.threshHiFreq > 0 && isfield(dPARAMS,'plotaxes')
        xtline = [dPARAMS.plotaxes.minRMS,dPARAMS.plotaxes.maxRMS];
        ytline = [p.threshHiFreq ,p.threshHiFreq];
    elseif p.threshHiFreq > 0
        xtline = [min(dHANDLES.pxmsp),max(dHANDLES.pxmsp)];
        ytline = [p.threshHiFreq ,p.threshHiFreq];
    end
    plot(dHANDLES.h53,xtline,ytline,'r')
end

if dPARAMS.ff2 % false in red
    plot(dHANDLES.h53,dPARAMS.pxmsp(dPARAMS.K2),dPARAMS.freq(dPARAMS.K2),...
        'r.','MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K2))
end
if dPARAMS.ff3 % ID'd in associated color
    
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        thisIDset = dPARAMS.spCodeSet == dPARAMS.specIDs(iC2);
        hPP = plot(dHANDLES.h53,dPARAMS.pxmsp(dPARAMS.K3(thisIDset)),...
            dPARAMS.freq(dPARAMS.K3(thisIDset)),'.',...
            'MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K3(thisIDset)));
        set(hPP,'Color',p.colorTab(dPARAMS.specIDs(iC2),:))
    end
end


hold(dHANDLES.h53, 'off')
if p.threshHiFreq > 0
    ylim(dHANDLES.h53,[p.fLow dPARAMS.ymax+10])
else
    ylim(dHANDLES.h53,[p.fLow dPARAMS.ymax])
end

xlabel(dHANDLES.h53,'dB RMS')
ylabel(dHANDLES.h53,'Peak Frequency (kHz)')
if isfield(dPARAMS,'plotaxes')
    xlim(dHANDLES.h53,[dPARAMS.plotaxes.minRMS,dPARAMS.plotaxes.maxRMS])
end

if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h53,'on')
    plot(dHANDLES.h53,dPARAMS.pxmsp(dPARAMS.yell),dPARAMS.freq(dPARAMS.yell),'ko',...
        'MarkerSize',p.sizeBack,...
        'LineWidth',2,'UserData',dPARAMS.clickTimes(dPARAMS.K2))
    hold(dHANDLES.h53,'off')
end