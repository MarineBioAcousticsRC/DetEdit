function highlight_brushed(ax,x,y,yellIdx)

% highlight_brushed
%
% Show yellow-brushed (highlighted) detections clearly:
%   1) redraw the highlighted points on top of everything else, in the
%      same style they normally have (unlabeled / false / ID'd, including
%      outline = prediction, fill = ground truth when plotCorrect_incorrect
%      is on), respecting the legend on/off toggles
%   2) circle each of them in black (a ring just outside the point, so the
%      points' own outline/fill colors stay fully visible)
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
ringSize = max(round(p.sizePoints*0.6),5); % size of the ID'd circles (plot_AR_rings)
circleSize = max(dotSize,ringSize) + 5; % black circle just outside the point
circleWidth = 1.5; % thickness of the black circle

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
        if isfield(dPARAMS,'showAR') && dPARAMS.showAR
            h = plot_AR_rings(ax,idX(inCode),idY(inCode),idT(inCode),...
                dPARAMS.arStatus(idLocs(inCode)),thisColor,p.sizePoints,...
                dPARAMS.ID_Toggle{thisCode});
            set(h,'Tag','yellFront','UserData',[]) % not brushable
        elseif dPARAMS.showGT
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

% black circle around each highlighted point that is currently shown
% (respecting the legend toggles). No 'UserData', so brushing ignores them.
isShown = (isUnlabeled & strcmp(dPARAMS.NoLabel_Toggle,'on')) | ...
    (isFD & strcmp(dPARAMS.FD_Toggle,'on'));
if any(isID)
    idRows = find(isID);
    for iR = 1:length(idRows)
        thisCode = dPARAMS.spCodeSet(idLoc(idRows(iR)));
        isShown(idRows(iR)) = strcmp(dPARAMS.ID_Toggle{thisCode},'on');
    end
end
xl = get(ax,'XLim'); yl = get(ax,'YLim');
if any(isShown)
    plot(ax,x(isShown),y(isShown),'o','MarkerSize',circleSize,...
        'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',circleWidth,...
        'LineStyle','none','Tag','yellFront');
end
set(ax,'XLim',xl,'YLim',yl) % keep the axis limits as they were

if ~holdState
    hold(ax,'off')
end
