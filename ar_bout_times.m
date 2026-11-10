function boutTimes = ar_bout_times(clickTimes, clickLevels, threshRL, entry)

% ar_bout_times
%
% Accept/reject review: the click times that define sessions (bouts).
% Only clicks that have an accepted or rejected label, and that pass the
% RL threshold, are used, so every session has something you reviewed in
% it. All clicks are still shown inside each session.
%
% mkLTSAsessions and detEdit both call this, so their session lists match.
% The labels come from the original accepted / rejected ID files (not the
% work-in-progress file), so editing labels doesn't change the sessions.
%
% Inputs:
%   clickTimes, clickLevels - TPWS click times and levels (tf already added)
%   threshRL - p.threshRL
%   entry - one element of ar_build_index (fields idAcc, idRej)
%
% Output:
%   boutTimes - sorted unique column of click times

labelTimes = [];
for f = {entry.idAcc, entry.idRej}
    if ~isempty(f{1}) && exist(f{1},'file') == 2
        S = load(f{1},'zID');
        if isfield(S,'zID') && ~isempty(S.zID)
            labelTimes = [labelTimes; double(S.zID(:,1))]; %#ok<AGROW>
        end
    end
end
clickTimes = clickTimes(:);
clickLevels = clickLevels(:);
keep = clickLevels >= threshRL & ismember(clickTimes,labelTimes);
boutTimes = unique(clickTimes(keep));
