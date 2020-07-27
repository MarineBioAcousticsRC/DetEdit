function figure51

global dPARAMS dHANDLES p

figure(dHANDLES.RMSvPPfig);clf
dHANDLES.h51 = gca;

if p.specploton && p.loadMSP
    plot(dHANDLES.h51,dPARAMS.transfRMSall,dPARAMS.xmppAll,'o','MarkerSize',...
        p.sizeBack,'MarkerEdgeColor',[.7,.7,.7],'UserData',dPARAMS.clickTimes)
    title(dHANDLES.h51,['Based on total of ',num2str(length(dPARAMS.xmppAll)),' clicks']);
    
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
dHANDLES.RL51 = plot(dHANDLES.h51,dPARAMS.transfRMS(dPARAMS.unlabeledIdx),...
    dPARAMS.xmpp(dPARAMS.unlabeledIdx),'.','MarkerSize',p.sizePoints,...
    'UserData',dPARAMS.t(dPARAMS.unlabeledIdx),'Visible',dPARAMS.NoLabel_Toggle);% true ones in blue


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


dHANDLES.RLFD51 = [];
if dPARAMS.ff2 % false in red
    dHANDLES.RLFD51 = plot(dHANDLES.h51,dPARAMS.transfRMS(dPARAMS.K2),...
        dPARAMS.xmpp(dPARAMS.K2),'r.','MarkerSize',p.sizePoints,...
        'UserData',dPARAMS.t(dPARAMS.K2),'Visible',dPARAMS.FD_Toggle);
end

dHANDLES.RLID51 = cell(size(p.colorTab,1),1);
if dPARAMS.ff3 % ID'd in associated color
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        iColor = dPARAMS.specIDs(iC2);
        dHANDLES.RLID51{iColor} =  plot(dHANDLES.h51,dPARAMS.transfRMS(dPARAMS.K3(dPARAMS.thisIDset{iC2})),...
            dPARAMS.xmpp(dPARAMS.K3(dPARAMS.thisIDset{iC2})),'.',...
            'MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.K3(dPARAMS.thisIDset{iC2})));
        set(dHANDLES.RLID51{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC2),:),...
            'Visible',dPARAMS.ID_Toggle{iColor})
        
    end
end
 
hold(dHANDLES.h51, 'off')

xlabel(dHANDLES.h51,'Transformed Received Level (dB_r_m_s re 1\muPa)')
ylabel(dHANDLES.h51,'Received Level (dB_p_p re 1\muPa)')
xlim(dHANDLES.h51,[p.rmsLow,p.rmsHi])
ylim(dHANDLES.h51,[p.rlLow,p.rlHi])


% if you have items brushed in yellow, highlight those on each plot
if ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h51,'on')
    plot(dHANDLES.h51,dPARAMS.transfRMS(dPARAMS.yell),dPARAMS.xmpp(dPARAMS.yell),'ko',...
        'MarkerSize',p.sizeBack,...
        'LineWidth',2,'UserData',dPARAMS.clickTimes(dPARAMS.K2))
    hold(dHANDLES.h51,'off')
end
