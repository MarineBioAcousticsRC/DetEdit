function [brushDate, brushColor, bFlag] = get_brushed(hFig)


brushDate = [];
brushColor = [];

% Find brushed axes
hBrushAll = get(hFig,'Children');
selLine = arrayfun(@(x) hBrushAll(x).Marker == '.', 1:length(hBrushAll), 'UniformOutput', false);
if iscell(selLine)
    lenAx = cell2mat(cellfun(@length,selLine,'UniformOutput',false));
    selLine(lenAx > 1) = {false};
    selLine = cell2mat(selLine);
end
brushedLine = find(selLine,1,'last');
hBrush = hBrushAll(brushedLine);

% Get indices of brushed data
brushID = find(get(hBrush, 'BrushData'));

bFlag = 0;

if ~isempty(brushID)
    
    bFlag = 1;
    
    % Get color info from brush
    h = get(brush);
    brushColor = round(h.Color*100)./100;
    
    % Get vector of dates associated with these points
    markerDates = hBrush.UserData;
    brushDate = markerDates(brushID);
end