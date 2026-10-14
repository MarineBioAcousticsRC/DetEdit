function highlight_brushed(ax,x,y,yellIdx)

% highlight_brushed
%
% Show yellow-brushed (highlighted) detections clearly:
%   1) redraw the highlighted points on top of everything else, in the
%      same style they normally have (unlabeled / false / ID'd, including
%      outline = prediction, fill = ground truth when plotCorrect_incorrect
%      is on), respecting the legend on/off toggles
%   2) draw one thick black box around all of them (so the points' own
%      outline/fill colors stay fully visible)
%
% The redrawn copies use the 'o' marker and a 'yellFront' tag, so brushing
% (get_brushed) ignores them.
%
% Inputs:
%   ax - axes to draw in
%   x, y - coordinates of the highlighted points
%   yellIdx - their indices into dPARAMS.t (same length as x and y)

global dPARAMS p

x = x(:); y = y(:); yellIdx = yellIdx(:);
if isempty(x)
    return
end
holdState = ishold(ax);
hold(ax,'on')

dotSize = max(round(p.sizePoints/2.5),3); % 'o' about the size of a '.' point
boxLineWidth = 2; % thickness of the box around the selection
boxPad = 0.015; % padding around the selection, as a fraction of the axis range

% sort points into false / ID'd / unlabeled, the same way boutMotion does
isFD = false(size(yellIdx));
if dPARAMS.ff2
    isFD = ismember(yellIdx,dPARAMS.K2);
end
isID = false(size(yellIdx));
idLoc = zeros(size(yellIdx));
if dPARAMS.ff3
    [isID,idLoc] = ismember(yellIdx,dPARAMS.K3);
    isID(isID) = dPARAMS.labelConfIdx(idLoc(isID)); % low confidence = unlabeled
    isID = isID & ~isFD;
end
isUnlabeled = ~isFD & ~isID;

% unlabeled
defaultColor = get(ax,'ColorOrder');
plot(ax,x(isUnlabeled),y(isUnlabeled),'o','MarkerSize',dotSize,...
    'MarkerEdgeColor',defaultColor(1,:),'MarkerFaceColor',defaultColor(1,:),...
    'LineStyle','none','Tag','yellFront','Visible',dPARAMS.NoLabel_Toggle);

% false
plot(ax,x(isFD),y(isFD),'o','MarkerSize',dotSize,'MarkerEdgeColor','r',...
    'MarkerFaceColor','r','LineStyle','none','Tag','yellFront',...
    'Visible',dPARAMS.FD_Toggle);

% ID'd, one group per label
if any(isID)
    codes = dPARAMS.spCodeSet(idLoc(isID));
    idX = x(isID); idY = y(isID); idLocs = idLoc(isID); idT = dPARAMS.t(yellIdx(isID));
    for thisCode = unique(codes(:))'
        inCode = codes == thisCode;
        thisColor = p.colorTab(thisCode,:);
        if dPARAMS.showGT
            h = plot_misclass_rings(ax,idX(inCode),idY(inCode),...
                idT(inCode),dPARAMS.gtCode(idLocs(inCode)),...
                thisColor,p.colorTab,p.sizePoints,dPARAMS.ID_Toggle{thisCode});
            set(h,'Tag','yellFront','UserData',[]) % not brushable
        else
            plot(ax,idX(inCode),idY(inCode),'o','MarkerSize',dotSize,...
                'MarkerEdgeColor',thisColor,'MarkerFaceColor',thisColor,...
                'LineStyle','none','Tag','yellFront',...
                'Visible',dPARAMS.ID_Toggle{thisCode});
        end
    end
end

% thick black box around all highlighted points (no markers, so brushing
% ignores it). Padding is relative to the current axis limits.
% Only points inside the visible plot area count, so an off-scale point
% can't stretch the box off the plot.
xl = get(ax,'XLim'); yl = get(ax,'YLim');
inView = x >= xl(1) & x <= xl(2) & y >= yl(1) & y <= yl(2);
if any(inView)
    xPad = boxPad*diff(xl); yPad = 2*boxPad*diff(yl);
    bx = [min(x(inView))-xPad, max(x(inView))+xPad];
    by = [min(y(inView))-yPad, max(y(inView))+yPad];
    plot(ax,bx([1,2,2,1,1]),by([1,1,2,2,1]),'k-','LineWidth',boxLineWidth,...
        'Marker','none','Tag','yellFront');
end
set(ax,'XLim',xl,'YLim',yl) % don't let the box change the axis limits

if ~holdState
    hold(ax,'off')
end
