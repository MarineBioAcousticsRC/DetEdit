function [brushDate, brushColor, bFlag, brushXY] = get_brushed(hFig)

% brushXY - [x y] plot coordinates of the brushed points (used to box the
%           exact ICI points brushed in the bottom LTSA panel)


brushDate = [];
brushColor = [];
brushXY = zeros(0,2);

% Find brushed axes
hBrushAll = get(hFig,'Children');
% brushable detection lines: '.' marker, or ID circles tagged 'IDring'
% (plotCorrect_incorrect mode). Objects without a Marker/BrushData property
% (images, text, ...) are skipped, which newer MATLAB versions need.
selLine = false(1,length(hBrushAll));
for iH = 1:length(hBrushAll)
    hThis = hBrushAll(iH);
    if ~isprop(hThis,'BrushData')
        continue
    end
    isDot = isprop(hThis,'Marker') && strcmp(get(hThis,'Marker'),'.');
    isRing = strcmp(get(hThis,'Tag'),'IDring');
    selLine(iH) = isDot || isRing;
end
brushedLine = find(selLine);%,1,'last');
hBrush = hBrushAll(brushedLine);

% Get indices of brushed data
brushID = [];
bFlag = 0;
for iB = 1:size(brushedLine,2)
    brushID = find(get(hBrush(iB), 'BrushData'));
    
    if ~isempty(brushID)
        bFlag = 1;
        
        % Get color info from brush
        try
            hb = brush(ancestor(hFig,'figure'));
            brushColor = round(get(hb,'Color')*100)./100;
        catch
            h = get(brush);
            brushColor = round(h.Color*100)./100;
        end
        
        % Get vector of dates associated with these points
        markerDates = hBrush(iB).UserData;
        brushDate = [brushDate;markerDates(brushID)];
        xd = get(hBrush(iB),'XData'); yd = get(hBrush(iB),'YData');
        xd = xd(:); yd = yd(:);
        brushXY = [brushXY; xd(brushID), yd(brushID)]; %#ok<AGROW>

    end
    
end