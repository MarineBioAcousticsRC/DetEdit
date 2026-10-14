function colorTab = fix_colorTab(colorTab,nNeeded)

% fix_colorTab
%
% Clean up the ID color table and make sure it has at least nNeeded rows.
%   - rows that are white or nearly white are replaced
%   - rows that exactly match a reserved brush color (red = false,
%     green = true, yellow = highlight) are replaced
%   - rows that duplicate an earlier row are replaced
%   - new rows are added if there are more labels than colors
% Replacement/new colors are picked to be as different as possible from
% every color already in use (and from white, black, the gray background
% points and the default blue used for unlabeled detections).
% Colors you set yourself are kept unless they break one of the rules above.
%
% Inputs:
%   colorTab - [N x 3] RGB table (0-1)
%   nNeeded - minimum number of rows required
%
% Output:
%   colorTab - [max(N,nNeeded) x 3] cleaned table, rounded to 2 decimals
%              (brushing matches colors at 2 decimals)

if nargin < 2 || isempty(nNeeded)
    nNeeded = 0;
end
colorTab = round(colorTab*100)/100;
nRows = max(size(colorTab,1),nNeeded);

reserved = [1,0,0; 0,1,0; 1,1,0]; % false, true, highlight brushes
% colors new picks should stay away from (not banned for user rows)
avoid = [1,1,1; 0,0,0; .7,.7,.7; 0,.45,.74; reserved];

% find rows that need replacing
bad = false(nRows,1);
bad(size(colorTab,1)+1:end) = true; % missing rows
for iRow = 1:size(colorTab,1)
    c = colorTab(iRow,:);
    isWhite = min(c) >= 0.85;
    isReserved = any(all(abs(reserved - c) < 0.005,2));
    isDup = iRow > 1 && any(all(abs(colorTab(1:iRow-1,:) - c) < 0.005,2));
    bad(iRow) = isWhite || isReserved || isDup;
end
if ~any(bad)
    return
end

% candidate pool: a spread of hues, saturations and brightnesses,
% skipping anything too light to see on a white background
[hh,ss,vv] = ndgrid(0:1/48:1-1/48,[.45,.7,1],[.45,.7,.9]);
pool = round(hsv2rgb([hh(:),ss(:),vv(:)])*100)/100;
lum = pool*[.2126;.7152;.0722];
pool = pool(min(pool,[],2) < 0.75 & lum < 0.85,:);
pool = unique(pool,'rows','stable');

% greedily pick the candidate farthest from everything in use
keep = colorTab(~bad(1:size(colorTab,1)),:);
inUse = [keep; avoid];
newTab = zeros(nRows,3);
newTab(1:size(colorTab,1),:) = colorTab;
for iRow = find(bad)'
    if isempty(pool)
        error('fix_colorTab: ran out of distinct colors (asked for %d labels).',nRows)
    end
    dists = inf(size(pool,1),1);
    for iUse = 1:size(inUse,1)
        dists = min(dists,sqrt(sum((pool - inUse(iUse,:)).^2,2)));
    end
    [~,best] = max(dists);
    newTab(iRow,:) = pool(best,:);
    inUse = [inUse; pool(best,:)]; %#ok<AGROW>
    pool(best,:) = [];
end
colorTab = newTab;
