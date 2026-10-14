function newRows = make_manual_ID_rows(times,label,nCols)

% make_manual_ID_rows
%
% Build zID rows for detections labeled by hand in detEdit (brushing,
% number keys, 'l', 'e'). A manual label is treated as confirmed, so:
%   col 1 = detection time
%   col 2 = label
%   col 3 = confidence -> 1
%   col 4 = ground truth -> label
%   col 5+ (if any) = NaN
%
% Inputs:
%   times - column of detection times
%   label - scalar label code
%   nCols - number of columns in zID (rows are built to match; 2 minimum)
%
% Output:
%   newRows - [length(times) x nCols] matrix to append to zID

times = times(:);
nCols = max(nCols,2);
nRows = length(times);

newRows = nan(nRows,nCols);
newRows(:,1) = times;
newRows(:,2) = label;
if nCols >= 3
    newRows(:,3) = 1;
end
if nCols >= 4
    newRows(:,4) = label;
end
