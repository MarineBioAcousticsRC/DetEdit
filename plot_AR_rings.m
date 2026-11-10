function h = plot_AR_rings(ax,x,y,t,status,labelColor,sizePoints,visibleState)

% plot_AR_rings
%
% Draw ID'd points as circles in accept/reject review mode
% (p.plotAcceptReject):
%   - filled center = color of the label (zID col 2)
%   - outline = same label color if accepted (so it looks like a solid dot),
%               red if rejected
%
% Same look and behavior as plot_misclass_rings: lines are tagged 'IDring'
% so get_brushed treats them like normal '.' detection points (UserData
% holds detection times), and each line gets a 'gtMatch' appdata value
% ('correct' = accepted, 'wrong' = rejected) so stack_ID_rings draws the
% rejected ones on top.
%
% Inputs:
%   ax - axes to draw in
%   x, y - point coordinates
%   t - detection times for each point (used by brushing)
%   status - 1 = accepted, 0 = rejected, for each point
%   labelColor - RGB color of the label
%   sizePoints - p.sizePoints
%   visibleState - 'on' or 'off'
%
% Output:
%   h - column of line handles (empty if nothing to plot)

ringWidth = 1.5; % outline thickness
ringSize = max(round(sizePoints*0.6),5); % circle diameter (points)
rejectColor = [1,0,0];

h = gobjects(0,1);
if isempty(x)
    return
end
x = x(:); y = y(:); t = t(:); status = status(:);
isRej = status == 0;

holdState = ishold(ax);
hold(ax,'on')
groups = {~isRej, isRej};
edgeColors = {labelColor, rejectColor};
matchNames = {'correct','wrong'};
for iG = 1:2
    sel = groups{iG};
    if ~any(sel)
        continue
    end
    h(end+1,1) = plot(ax,x(sel),y(sel),'o',...
        'MarkerSize',ringSize,'LineWidth',ringWidth,...
        'MarkerEdgeColor',edgeColors{iG},'MarkerFaceColor',labelColor,...
        'LineStyle','none','UserData',t(sel),'Tag','IDring',...
        'Visible',visibleState); %#ok<AGROW>
    setappdata(h(end),'gtMatch',matchNames{iG})
end
if ~holdState
    hold(ax,'off')
end
