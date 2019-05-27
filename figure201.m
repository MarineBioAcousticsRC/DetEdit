function figure201

global dPARAMS dHANDLES fNameList p

warning('off') % why turn off warning? why not fix problem?
figure(dHANDLES.LTSAfig);clf
dHANDLES.LTSAsubs = subplot_layout; % Top panel, Figure 201: Received Level

% set(dHANDLES.hbLTSA,'ActionPostCallback',@brushAction)

plot(dHANDLES.LTSAsubs(1),dPARAMS.t,dPARAMS.RL,'.','MarkerSize',p.sizePoints,...
    'UserData',dPARAMS.t)
hold(dHANDLES.LTSAsubs(1),'on')
if dPARAMS.ff2 % plot False detections in red
    plot(dHANDLES.LTSAsubs(1),dPARAMS.tfd,dPARAMS.rlFD,'r.','MarkerSize',p.sizePoints,'UserData',dPARAMS.tfd)
    % disp([' false det plotted:',num2str(length(tfd))])
end
if dPARAMS.ff3 % plot ID'd detections in associated color
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        thisIDset = dPARAMS.spCodeSet ==dPARAMS.specIDs(iC2);
        hRLID = plot(dHANDLES.LTSAsubs(1),dPARAMS.tID(thisIDset),dPARAMS.rlID(thisIDset),'.',...
            'MarkerSize',p.sizePoints,'UserData',dPARAMS.tID(thisIDset));
        set(hRLID,'Color',p.colorTab(dPARAMS.specIDs(iC2),:))
    end
end
hold(dHANDLES.LTSAsubs(1),'off')
axis(dHANDLES.LTSAsubs(1),[dPARAMS.PT(1) dPARAMS.PT(end) p.rlLow p.rlHi])
datetick(dHANDLES.LTSAsubs(1),'x',15,'keeplimits')
grid(dHANDLES.LTSAsubs(1),'on')
tstr(1) = {fNameList.TPWS};
tstr(2) = {['Session: ',num2str(dPARAMS.k),'/',num2str(dPARAMS.nb),' Start Time ',...
    datestr(dPARAMS.sb(dPARAMS.k)),' Detect = ',num2str(dPARAMS.nd)]};
title(dHANDLES.LTSAsubs(1),tstr);
ylabel(dHANDLES.LTSAsubs(1),'RL [dB re 1\muPa]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% middle panel LTSA
c = (p.ltsaContrast/100) .* dPARAMS.pwr1 + p.ltsaBright;
image(dPARAMS.PT,dPARAMS.f/1000,c,'parent',dHANDLES.LTSAsubs(2))
set(dHANDLES.LTSAsubs(2),'yDir','normal')
axis(dHANDLES.LTSAsubs(2),[dPARAMS.PT(1) dPARAMS.PT(end) p.ltsaLims])%v2(4)
ylabel(dHANDLES.LTSAsubs(2),'Frequency (kHz)')
datetick(dHANDLES.LTSAsubs(2),'keeplimits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom panel, Figure 201: Inter-Detection Interval
% make two copies of dt points for brush
tdt2 = [];
dt2 = [];
ldt = length(dPARAMS.dt);
if ldt > 0
    tdt2 = reshape([dPARAMS.t(1:ldt),dPARAMS.t((1:ldt)+1)]',2*ldt,1);
    dt2 = reshape([dPARAMS.dt,dPARAMS.dt]',2*ldt,1);
    
    hICI = plot(dHANDLES.LTSAsubs(3),tdt2,dt2,'.');
    set(hICI,'MarkerSize',p.sizePoints,'MarkerFaceColor','b','LineStyle','none','UserData',tdt2)
    
    % Do setup for 1st axes
    axis(dHANDLES.LTSAsubs(3),[dPARAMS.PT(1) dPARAMS.PT(end) 0 p.dtHi])
    datetick(dHANDLES.LTSAsubs(3),'x',15,'keeplimits')
    Ytick = 0:p.dtHi/10:p.dtHi;
    set(dHANDLES.LTSAsubs(3),'YTick',Ytick)
    datetick(dHANDLES.LTSAsubs(3),'x',15,'keeplimits')
    grid(dHANDLES.LTSAsubs(3),'on')
    ylabel(dHANDLES.LTSAsubs(3),'Time between detections [s]')
    
    %%% plot FD, ID
    hold(dHANDLES.LTSAsubs(3),'on')
    if dPARAMS.ff2
        plot(dHANDLES.LTSAsubs(3), dPARAMS.tfd(2:end), dPARAMS.dtFD,'.r',...
            'MarkerSize',p.sizePoints,'UserData', dPARAMS.tfd(2:end))
        % no need to double FD since only the blue points are brush captured
    end
    if dPARAMS.ff3 % plot ID'd in associated color
        for iC3 = 1:length(dPARAMS.specIDs) % set colors
            thisIDset = dPARAMS.spCodeSet == dPARAMS.specIDs(iC3);
            hICI = plot(dHANDLES.LTSAsubs(3),dPARAMS.tID(thisIDset(2:end)),...
                dPARAMS.dtID(thisIDset(2:end)),'.','MarkerSize',p.sizePoints,...
                'UserData',dPARAMS.tID(thisIDset));
            set(hICI,'Color',p.colorTab(dPARAMS.specIDs(iC3),:))
        end
    end
    hold(dHANDLES.LTSAsubs(3),'off')
else
    plot(0,0);
end

if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    hold(dHANDLES.LTSAsubs(1),'on')
    plot(dHANDLES.LTSAsubs(1),dPARAMS.t(dPARAMS.yell),dPARAMS.RL(dPARAMS.yell),...
        'ko','MarkerSize',p.sizeBack,'UserData',dPARAMS.t(dPARAMS.yell));
    hold(dHANDLES.LTSAsubs(1),'off');
    
    
    % for diffs, yell can't exceed length dt, which could happen if you
    % grabbed the last point in the vector, so:
    yellDT = dPARAMS.yell(dPARAMS.yell<length(dPARAMS.dt));
    hold(dHANDLES.LTSAsubs(3),'on')
    plot(dHANDLES.LTSAsubs(3),dPARAMS.t(yellDT),dPARAMS.dt(yellDT),'ko','MarkerSize',...
        p.sizeBack,'UserData',dPARAMS.t(dPARAMS.yell));
    hold(dHANDLES.LTSAsubs(3),'off')
    
end