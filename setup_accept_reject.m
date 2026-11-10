function fNameList = setup_accept_reject(fNameList, p)

% setup_accept_reject
%
% Accept/reject review mode (p.plotAcceptReject = true).
%
% Reads the reviewed labels for this disk from its accepted and rejected
% ID files (p.arEntry, from ar_build_index: matched by deployment + disk)
%
% and builds one working ID file that detEdit edits as usual:
%   <outDir>\checked\inProgress\<ID file>
% with zID columns
%   1 = time, 2 = label, 3 = confidence, 4 = status (1 = accepted, 0 = rejected)
%
% If a working file from an earlier run exists, you're asked whether to
% resume it or start over from the accepted / rejected files.
%
% The final file (accepted only, 3 columns) is written by save_checked_ID
% to <outDir>\checked\<ID file>.
%
% Output fNameList gets these fields:
%   ID        - working file (what detEdit loads/saves while you edit)
%   IDacc     - accepted file ('' if none)
%   IDrej     - rejected file ('' if none)
%   IDchecked - final output file

e = p.arEntry;
idfn = e.idName;
fNameList.IDacc = e.idAcc;
fNameList.IDrej = e.idRej;
checkedDir = fullfile(p.outDir,'checked');
workDir = fullfile(checkedDir,'inProgress');
fNameList.IDchecked = fullfile(checkedDir,idfn);
fNameList.ID = fullfile(workDir,idfn);

hasAcc = ~isempty(fNameList.IDacc) && exist(fNameList.IDacc,'file') == 2;
hasRej = ~isempty(fNameList.IDrej) && exist(fNameList.IDrej,'file') == 2;

fprintf('\nAccept/reject review mode: %s\n',e.label);
fprintf('  accepted: %s%s\n',fNameList.IDacc,missing_str(hasAcc));
fprintf('  rejected: %s%s\n',fNameList.IDrej,missing_str(hasRej));
fprintf('  checked : %s\n',fNameList.IDchecked);

if ~exist(workDir,'dir')
    mkdir(workDir);
end

% resume an earlier session?
if exist(fNameList.ID,'file') == 2
    fprintf('  Found work in progress: %s\n',fNameList.ID);
    if exist(fNameList.IDchecked,'file') == 2
        fprintf('  (this disk was already saved to checked\\ before)\n');
    end
    answer = '';
    while ~any(strcmpi(answer,{'r','s'}))
        answer = strtrim(input('  Resume it (r) or start over from accepted/rejected (s)? ','s'));
    end
    if strcmpi(answer,'r')
        fprintf('  Resuming work in progress.\n\n');
        return
    end
end

if ~hasAcc && ~hasRej
    error('detEdit:acceptReject',['No accepted or rejected ID file found for this disk:\n',...
        '  %s\n  %s'],fNameList.IDacc,fNameList.IDrej);
end

% read both files
[zAcc,Sacc] = read_ID(fNameList.IDacc,hasAcc);
[zRej,Srej] = read_ID(fNameList.IDrej,hasRej);

% species list: from the accepted file, else the rejected one
mySpID = [];
if isfield(Sacc,'mySpID')
    mySpID = Sacc.mySpID;
elseif isfield(Srej,'mySpID')
    mySpID = Srej.mySpID;
end
if hasAcc && hasRej && isfield(Sacc,'mySpID') && isfield(Srej,'mySpID') && ...
        ~isequal(Sacc.mySpID,Srej.mySpID)
    warning('detEdit:acceptReject',['mySpID differs between the accepted and ',...
        'rejected files; using the accepted file''s list.'])
end

% accepted first, so a time that is (unexpectedly) in both stays accepted
zID = [zAcc, ones(size(zAcc,1),1); zRej, zeros(size(zRej,1),1)];
if ~isempty(zID)
    [~,firstRow] = unique(zID(:,1),'first');
    nDup = size(zID,1) - length(firstRow);
    if nDup > 0
        fprintf('  [WARN] %d time(s) in both accepted and rejected: kept as accepted.\n',nDup);
    end
    zID = zID(sort(firstRow),:);
    zID = sortrows(zID,1);
else
    zID = zeros(0,4);
end
fprintf('  Loaded %d accepted and %d rejected labels.\n\n',...
    sum(zID(:,4) == 1),sum(zID(:,4) == 0));

save(fNameList.ID,'zID','mySpID');

end

% ----------
function [z,S] = read_ID(fn,exists)
% zID of one file as [time, label, confidence] (confidence NaN if missing)
z = zeros(0,3);
S = struct();
if ~exists
    return
end
S = load(fn);
if ~isfield(S,'zID') || isempty(S.zID)
    return
end
z = double(S.zID);
if size(z,2) < 3
    z = [z, nan(size(z,1),3-size(z,2))];
end
z = z(:,1:3);
end

% ----------
function s = missing_str(exists)
if exists
    s = '';
else
    s = '   (not found)';
end
end
