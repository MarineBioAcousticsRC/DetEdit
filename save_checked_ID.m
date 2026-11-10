function save_checked_ID(zID, fNameList)

% save_checked_ID
%
% Accept/reject review mode (p.plotAcceptReject): write the final ID file
% for this disk to <outDir>\checked\<ID file>.
%
% Only ACCEPTED labels are kept (zID col 4 == 1; a row without a status
% counts as accepted). Rejected labels you didn't flip are left out
% entirely. The file has the usual 3 zID columns: time, label, confidence.
%
% All other variables (mySpID etc.) are copied from the accepted file (or
% the rejected file if there is no accepted one), like the review script
% does. The work-in-progress file in checked\inProgress\ is kept, so the
% disk can be reopened later with its accepted/rejected state intact.
%
% Inputs:
%   zID       - working ID matrix [time, label, confidence, status]
%   fNameList - struct with fields IDchecked, IDacc, IDrej (setup_accept_reject)

if isempty(zID)
    zID = zeros(0,4);
end
if size(zID,2) < 3
    zID(:,end+1:3) = NaN;
end
if size(zID,2) >= 4
    isAcc = zID(:,4) ~= 0; % 1 or NaN = accepted, 0 = rejected
else
    isAcc = true(size(zID,1),1);
end

% everything else in the file comes from the original accepted / rejected file
if exist(fNameList.IDacc,'file') == 2
    S = load(fNameList.IDacc);
elseif exist(fNameList.IDrej,'file') == 2
    S = load(fNameList.IDrej);
else
    S = struct();
end
S.zID = sortrows(zID(isAcc,1:3),1);
if ~isfield(S,'mySpID')
    W = load(fNameList.ID,'mySpID'); % working file
    if isfield(W,'mySpID')
        S.mySpID = W.mySpID;
    end
end

checkedDir = fileparts(fNameList.IDchecked);
if ~exist(checkedDir,'dir')
    mkdir(checkedDir);
end
save(fNameList.IDchecked,'-struct','S','-v7.3');

disp(' ')
fprintf('Saved checked ID file: %s\n',fNameList.IDchecked);
fprintf('  %d accepted labels saved, %d rejected labels left out\n',...
    sum(isAcc),sum(~isAcc));
