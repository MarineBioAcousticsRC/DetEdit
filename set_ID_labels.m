function zID = set_ID_labels(zID,times,label)

% set_ID_labels
%
% Apply a hand-picked label to detections, respecting the "Edit:" choice in
% the ID legend (dPARAMS.editTarget):
%
%   'both'  (default) - model output (col 2) AND ground truth (col 4) are set
%                       to the label, confidence (col 3) = 1
%   'model'           - only the model output (col 2) is set, confidence = 1;
%                       an existing ground truth is kept, new rows get none
%   'gt'              - only the ground truth (col 4) is set; the model
%                       output and confidence are kept. Detections with no
%                       model label yet get the label in both columns
%                       (a row can't exist without a model label).
%
% Inputs:
%   zID   - ID matrix [time, label, (confidence), (ground truth), ...]
%   times - detection times to label
%   label - label code
%
% Output:
%   zID   - updated matrix (widened to 4 columns for 'model'/'gt' editing)

global dPARAMS

editTarget = 'both';
if isfield(dPARAMS,'editTarget') && ~isempty(dPARAMS.editTarget)
    editTarget = dPARAMS.editTarget;
end

times = times(:);
if isempty(times)
    return
end

% make sure there's a ground truth column when editing one side only
nCols = size(zID,2);
if isempty(zID)
    nCols = max(nCols,2);
end
if ~strcmp(editTarget,'both') && nCols < 4
    zID = [zID, nan(size(zID,1),4-nCols)];
    nCols = 4;
end

% detections that already have a row
if isempty(zID)
    isOld = false(size(times)); loc = zeros(size(times));
else
    [isOld,loc] = ismember(times,zID(:,1));
end
rows = loc(isOld);
switch editTarget
    case 'model'
        zID(rows,2) = label;
        if nCols >= 3, zID(rows,3) = 1; end
    case 'gt'
        zID(rows,4) = label;
    otherwise % both
        zID(rows,2) = label;
        if nCols >= 3, zID(rows,3) = 1; end
        if nCols >= 4, zID(rows,4) = label; end
end

% detections without a row yet
newTimes = times(~isOld);
if ~isempty(newTimes)
    newRows = make_manual_ID_rows(newTimes,label,nCols);
    if strcmp(editTarget,'model') && nCols >= 4
        newRows(:,4) = NaN; % model edit only: no ground truth known
    elseif strcmp(editTarget,'gt')
        fprintf(['  %d detection(s) had no model label; gave them label %d ',...
            'as both model output and ground truth\n'],length(newTimes),label);
    end
    zID = [zID; newRows];
end
