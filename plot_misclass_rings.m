function h = plot_misclass_rings(ax,x,y,t,gtCode,predColor,colorTab,sizePoints,visibleState)

% plot_misclass_rings
%
% Draw ID'd points as circles when p.plotCorrect_incorrect is on:
%   - thick outline = color of the model's predicted label (zID col 2)
%   - filled center = color of the ground truth label (zID col 4)
%   - center left transparent when there is no ground truth (0 or NaN)
% So a correct label looks like a solid dot, and a wrong one shows two
% colors.
%
% One line is drawn per ground-truth class so each can have its own fill.
% Lines are tagged 'IDring' so get_brushed treats them like the normal '.'
% detection points (UserData holds detection times).
%
% Inputs:
%   ax - axes to draw in
%   x, y - point coordinates
%   t - detection times for each point (used by brushing)
%   gtCode - ground truth code for each point (0/NaN = none)
%   predColor - RGB color of the predicted label
%   colorTab - p.colorTab, used to look up ground truth colors
%   sizePoints - p.sizePoints
%   visibleState - 'on' or 'off'
%
% Output:
%   h - column of line handles (empty if nothing to plot)

ringWidth = 1.5; % outline thickness
ringSize = max(round(sizePoints*0.6),5); % circle diameter (points)

h = gobjects(0,1);
if isempty(x)
    return
end
x = x(:); y = y(:); t = t(:); gtCode = gtCode(:);

% anything without a usable ground truth color gets a transparent center
noGT = isnan(gtCode) | gtCode <= 0 | gtCode > size(colorTab,1) | ...
    gtCode ~= round(gtCode);
gtCode(noGT) = 0;

holdState = ishold(ax);
hold(ax,'on')
gtGroups = unique(gtCode);
for iG = 1:length(gtGroups)
    thisGroup = gtCode == gtGroups(iG);
    if gtGroups(iG) == 0
        faceColor = 'none';
    else
        faceColor = colorTab(gtGroups(iG),:);
    end
    h(end+1,1) = plot(ax,x(thisGroup),y(thisGroup),'o',...
        'MarkerSize',ringSize,'LineWidth',ringWidth,...
        'MarkerEdgeColor',predColor,'MarkerFaceColor',faceColor,...
        'LineStyle','none','UserData',t(thisGroup),'Tag','IDring',...
        'Visible',visibleState); %#ok<AGROW>
    % remember how this group compares to ground truth, for draw order
    % (see stack_ID_rings): correct / none / wrong
    if gtGroups(iG) == 0
        setappdata(h(end),'gtMatch','none')
    elseif isequal(faceColor,predColor)
        setappdata(h(end),'gtMatch','correct')
    else
        setappdata(h(end),'gtMatch','wrong')
    end
end
if ~holdState
    hold(ax,'off')
end
