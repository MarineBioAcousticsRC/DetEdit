function figure51

global dPARAMS dHANDLES p

figure(dHANDLES.RMSvPPfig);clf
dHANDLES.h51 = gca;

if p.specploton && p.loadMSP
    plot(dHANDLES.h51,dPARAMS.pxmspAll,dPARAMS.xmppAll,'o','MarkerSize',...
        p.sizeBack,'MarkerEdgeColor',[.7,.7,.7],'UserData',dPARAMS.clickTimes)
    title(dHANDLES.h51,['Based on ',num2str(length(dPARAMS.xmppAll)),' clicks']);
    
    if (p.threshRMS > 0)
        %%% why is x/ytline used for different variables?
        if p.threshPP > 0
            xtline = [p.threshRMS,p.threshRMS];
            ytline = [ min(dPARAMS.xmppAll),p.threshPP];
        else
            xtline = [p.threshRMS,p.threshRMS];
            ytline = [ min(dPARAMS.xmppAll),max(dPARAMS.xmppAll)];
        end
        hold(dHANDLES.h51,'on');
        plot(dHANDLES.h51,xtline,ytline,'r')
        hold(dHANDLES.h51,'off');
        
    end
end

% Plot  PP versus RMS Plot for this session
hold(dHANDLES.h51, 'on')
dHANDLES.h51.ColorOrderIndex = 1;
plot(dHANDLES.h51,dPARAMS.pxmsp,dPARAMS.xmpp,'.','MarkerSize',p.sizePoints,'UserData',dPARAMS.t)% true ones in blue

if p.threshPP > 0 && isfield(dPARAMS,'plotaxes') % why doesn't plot axes just either exist or not??? messy.
    xtline = [p.threshRMS,p.threshRMS];
    ytline = [dPARAMS.plotaxes.minPP,p.threshPP];
elseif p.threshPP > 0
    xtline = [p.threshRMS,p.threshRMS];
    ytline = [ min(dPARAMS.xmpp),p.threshPP]; % should default to same axes as PP timeseries
else
    xtline = [p.threshRMS,p.threshRMS];
    ytline = [min(dPARAMS.xmpp),max(dPARAMS.xmpp)];
end
plot(dHANDLES.h51,xtline,ytline,'r')


if dPARAMS.ff2 % false in red
    plot(dHANDLES.h51,dPARAMS.pxmsp(dPARAMS.K2),dPARAMS.xmpp(dPARAMS.K2),...
        'r.','MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K2))
end
if dPARAMS.ff3 % ID'd in associated color
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        thisIDset = dPARAMS.spCodeSet ==dPARAMS.specIDs(iC2);
        hPP = plot(dHANDLES.h51,dPARAMS.pxmsp(dPARAMS.K3(thisIDset)),...
            dPARAMS.xmpp(dPARAMS.K3(thisIDset)),'.',...
            'MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K3(thisIDset)));
        set(hPP,'Color',p.colorTab(dPARAMS.specIDs(iC2),:))
    end
end

hold(dHANDLES.h51, 'off')

xlabel(dHANDLES.h51,'dB RMS')
ylabel(dHANDLES.h51,'dB Peak-to-peak')
if isfield(dPARAMS,'plotaxes')
    xlim(dHANDLES.h51,[dPARAMS.plotaxes.minRMS,dPARAMS.plotaxes.maxRMS])
    ylim(dHANDLES.h51,[dPARAMS.plotaxes.minPP,dPARAMS.plotaxes.maxPP])
end


% if you have items brushed in yellow, highlight those on each plot
if ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h51,'on')
    plot(dHANDLES.h51,dPARAMS.pxmsp(dPARAMS.yell),dPARAMS.xmpp(dPARAMS.yell),'ko',...
        'MarkerSize',p.sizeBack,...
        'LineWidth',2,'UserData',dPARAMS.clickTimes(dPARAMS.K2))
    hold(dHANDLES.h51,'off')
end
