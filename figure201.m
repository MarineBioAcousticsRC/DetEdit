function figure201

global dPARAMS dHANDLES fNameList p

warning('off') % why turn off warning? why not fix problem?
figure(dHANDLES.LTSAfig);clf
dHANDLES.LTSAsubs = subplot_layout; % Top panel, Figure 201: Received Level

% set(dHANDLES.hbLTSA,'ActionPostCallback',@brushAction)

% only plot unlabeled stuff.
if ~isempty(dPARAMS.unlabeledIdx)
    dHANDLES.RL201 = plot(dHANDLES.LTSAsubs(1),dPARAMS.t(dPARAMS.unlabeledIdx),...
        dPARAMS.RL(dPARAMS.unlabeledIdx),...
        '.','MarkerSize',p.sizePoints,'UserData',dPARAMS.t(dPARAMS.unlabeledIdx),...
        'Visible',dPARAMS.NoLabel_Toggle);
end

dHANDLES.RLFD201 = [];
hold(dHANDLES.LTSAsubs(1),'on')
if dPARAMS.ff2 % plot False detections in red
    dHANDLES.RLFD201 = plot(dHANDLES.LTSAsubs(1),dPARAMS.tfd,dPARAMS.rlFD,...
        'r.','MarkerSize',p.sizePoints,'UserData',dPARAMS.tfd);
     set(dHANDLES.RLFD201,'Visible',dPARAMS.FD_Toggle)
    % disp([' false det plotted:',num2str(length(tfd))])
end

dHANDLES.RLID201 = cell(size(p.colorTab,1),1);
if dPARAMS.ff3 % plot ID'd detections in associated color
    for iC2 = 1:length(dPARAMS.specIDs) % set colors
        iColor = dPARAMS.specIDs(iC2);
        thisSet = dPARAMS.thisIDset{iC2};
        if dPARAMS.showGT % outline = prediction, fill = ground truth
            dHANDLES.RLID201{iColor} = plot_misclass_rings(dHANDLES.LTSAsubs(1),...
                dPARAMS.tID(thisSet),dPARAMS.rlID(thisSet),dPARAMS.tID(thisSet),...
                dPARAMS.gtCode(thisSet),p.colorTab(iColor,:),p.colorTab,...
                p.sizePoints,dPARAMS.ID_Toggle{iColor});
        else
            dHANDLES.RLID201{iColor} = plot(dHANDLES.LTSAsubs(1),dPARAMS.tID(thisSet),...
                dPARAMS.rlID(thisSet),'.',...
                'MarkerSize',p.sizePoints,'UserData',dPARAMS.tID(thisSet));
            set(dHANDLES.RLID201{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC2),:),...
                'Visible',dPARAMS.ID_Toggle{iColor})
        end
    end
end
if dPARAMS.ff3 && dPARAMS.showGT
    stack_ID_rings(dHANDLES.LTSAsubs(1)) % wrong on top, correct on bottom
end
hold(dHANDLES.LTSAsubs(1),'off')
axis(dHANDLES.LTSAsubs(1),[dPARAMS.PT(1) dPARAMS.PT(end) p.rlLow p.rlHi])
datetick(dHANDLES.LTSAsubs(1),'x',15,'keeplimits')
grid(dHANDLES.LTSAsubs(1),'on')
tstr(1) = {fNameList.TPWS};
tstr(2) = {['Session: ',num2str(dPARAMS.k),'/',num2str(dPARAMS.nb),' Start Time ',...
    datestr(dPARAMS.sb(dPARAMS.k)),' Detect = ',num2str(dPARAMS.nd)]};
if isfield(dPARAMS,'editTarget') && strcmp(dPARAMS.editTarget,'model')
    tstr{2} = [tstr{2},'   |   Labels edit: MODEL ONLY'];
elseif isfield(dPARAMS,'editTarget') && strcmp(dPARAMS.editTarget,'gt')
    tstr{2} = [tstr{2},'   |   Labels edit: GT ONLY'];
end
title(dHANDLES.LTSAsubs(1),tstr);
ylabel(dHANDLES.LTSAsubs(1),'Received Level (dB_p_p re 1\muPa)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% middle panel LTSA
c = (p.ltsaContrast/100) .* dPARAMS.pwr1 + p.ltsaBright;
image(dPARAMS.PT,dPARAMS.f/1000,c,'parent',dHANDLES.LTSAsubs(2))
set(dHANDLES.LTSAsubs(2),'yDir','normal')
axis(dHANDLES.LTSAsubs(2),[dPARAMS.PT(1) dPARAMS.PT(end) p.ltsaLims])%v2(4)
ylabel(dHANDLES.LTSAsubs(2),'Frequency (kHz)')
datetick(dHANDLES.LTSAsubs(2),'x',15,'keeplimits')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom panel, Figure 201: Inter-Detection Interval
% make two copies of dt points for brush
tdt2 = [];
dt2 = [];
ldt = length(dPARAMS.dt);
if ldt > 0
    % only plot unlabeled stuff, make duplicates
    if ~isempty(dPARAMS.unlabeledIdx)
        nUnlabeled = length(dPARAMS.unlabeledIdx)-1;
        tdt2 = reshape([dPARAMS.t(dPARAMS.unlabeledIdx(1:end-1)),...
            dPARAMS.t(dPARAMS.unlabeledIdx(2:end))]',2*nUnlabeled,1);
        dt2 = reshape([dPARAMS.dtUnlabeled,dPARAMS.dtUnlabeled]',2*(nUnlabeled),1);
    end
        
    
    dHANDLES.ICI201 = plot(dHANDLES.LTSAsubs(3),tdt2,dt2,'.');
    set(dHANDLES.ICI201,'MarkerSize',p.sizePoints,'MarkerFaceColor','b',...
        'LineStyle','none','UserData',tdt2,'Visible',dPARAMS.NoLabel_Toggle)
    
    % Do setup for 1st axes
    axis(dHANDLES.LTSAsubs(3),[dPARAMS.PT(1) dPARAMS.PT(end) 0 p.dtHi])
    datetick(dHANDLES.LTSAsubs(3),'x',15,'keeplimits')
    Ytick = 0:p.dtHi/10:p.dtHi;
    set(dHANDLES.LTSAsubs(3),'YTick',Ytick,'YTickLabelMode','auto')
    datetick(dHANDLES.LTSAsubs(3),'x',15,'keeplimits')
    grid(dHANDLES.LTSAsubs(3),'on')
    ylabel(dHANDLES.LTSAsubs(3),'Time between detections (s)')
    xlabel(dHANDLES.LTSAsubs(3),'Time (GMT)')
    
    %%% plot FD, ID
    hold(dHANDLES.LTSAsubs(3),'on')
    dHANDLES.ICIFD201 = [];
    if dPARAMS.ff2
        %duplicate again
        fdt2 = [dPARAMS.tfd(1:end-1);dPARAMS.tfd(2:end)];
        fddt2 = [dPARAMS.dtFD;dPARAMS.dtFD];
        dHANDLES.ICIFD201 = plot(dHANDLES.LTSAsubs(3),fdt2,fddt2,'.r',...
            'MarkerSize',p.sizePoints,'UserData',fdt2,'Visible',dPARAMS.FD_Toggle);
        % no need to double FD since only the blue points are brush captured
    end
    
    dHANDLES.ICIID201 = cell(size(p.colorTab,1),1);
    if dPARAMS.ff3 % plot ID'd in associated color
        for iC3 = 1:length(dPARAMS.specIDs) % set colors
            iColor = dPARAMS.specIDs(iC3);
            % limit to only things that meet the min confidence cutoff.
            thisIDsetTemp = find(dPARAMS.spCodeSet == dPARAMS.specIDs(iC3));
            thisIDset = thisIDsetTemp(dPARAMS.labelConfIdx(thisIDsetTemp));
            idt2 = [dPARAMS.tID(thisIDset(1:end-1));...
                dPARAMS.tID(thisIDset(2:end))];
            iddt2 = [dPARAMS.dtID(thisIDset(1:end-1));...
                dPARAMS.dtID(thisIDset(1:end-1))];
            
            if dPARAMS.showGT % outline = prediction, fill = ground truth
                igt2 = [dPARAMS.gtCode(thisIDset(1:end-1));...
                    dPARAMS.gtCode(thisIDset(2:end))];
                dHANDLES.ICIID201{iColor} = plot_misclass_rings(dHANDLES.LTSAsubs(3),...
                    idt2,iddt2,idt2,igt2,p.colorTab(iColor,:),p.colorTab,...
                    p.sizePoints,dPARAMS.ID_Toggle{iColor});
            else
                dHANDLES.ICIID201{iColor} = plot(dHANDLES.LTSAsubs(3),idt2,...
                     iddt2,'.','MarkerSize',p.sizePoints,'UserData',idt2);
                set(dHANDLES.ICIID201{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC3),:),...
                    'Visible',dPARAMS.ID_Toggle{iColor})
            end
        end
    end
    if dPARAMS.ff3 && dPARAMS.showGT
        stack_ID_rings(dHANDLES.LTSAsubs(3)) % wrong on top, correct on bottom
    end
    hold(dHANDLES.LTSAsubs(3),'off')
else
    plot(0,0);
end

if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    % top panel: bring highlighted points to the front and circle them
    highlight_brushed(dHANDLES.LTSAsubs(1),dPARAMS.t(dPARAMS.yell),...
        dPARAMS.RL(dPARAMS.yell),dPARAMS.yell);

    % LTSA: thin box spanning the highlighted points, +/- 2 minutes
    tPad = 2/(24*60); % 2 minutes in days
    tBox = [min(dPARAMS.t(dPARAMS.yell))-tPad, max(dPARAMS.t(dPARAMS.yell))+tPad];
    fLim = get(dHANDLES.LTSAsubs(2),'YLim');
    fBox = fLim + [1,-1]*0.02*diff(fLim); % inset slightly so top/bottom edges show
    hold(dHANDLES.LTSAsubs(2),'on')
    plot(dHANDLES.LTSAsubs(2),tBox([1,2,2,1,1]),fBox([1,1,2,2,1]),'-','Color',[1,1,1],...
        'LineWidth',1,'Tag','yellFront');
    hold(dHANDLES.LTSAsubs(2),'off')

    % bottom panel: if the points were brushed in this panel, box exactly
    % those. Otherwise mark the intervals (ICIs) *between* selected clicks.
    % Brushing one ICI point selects the two clicks on either end of it;
    % click j's own interval dt(j) runs to click j+1, so only keep j when
    % j+1 is selected too. Otherwise the second click's interval to the
    % next (unselected) click would be marked, which can be far off-scale.
    % If no two selected clicks are neighbours (e.g. single clicks picked
    % in the top panel), fall back to each click's interval to the next.
    if isfield(dPARAMS,'yellICIxy') && ~isempty(dPARAMS.yellICIxy)
        % brushed directly in this panel: use exactly the points brushed
        xy = unique(dPARAMS.yellICIxy,'rows');
        [onT,iT] = ismember(xy(:,1),dPARAMS.t);
        highlight_brushed(dHANDLES.LTSAsubs(3),xy(onT,1),xy(onT,2),iT(onT));
    else
        yellSel = dPARAMS.yell(:);
        nDT = length(dPARAMS.dt);
        yellDT = yellSel(ismember(yellSel+1,yellSel) & yellSel <= nDT);
        if isempty(yellDT)
            yellDT = yellSel(yellSel <= nDT);
        end
        highlight_brushed(dHANDLES.LTSAsubs(3),dPARAMS.t(yellDT),...
            dPARAMS.dt(yellDT),yellDT);
    end

end