function files = ar_build_index(p)

% ar_build_index
%
% Accept/reject review: find every disk that has an accepted and/or a
% rejected ID file, and match it to its TPWS file and LTSA files by
% deployment + disk (see ar_parse_name), not by full file name.
%
%   ID  <-> TPWS : deployment + disk + letter   (NFC_A_02 + disk01a)
%   TPWS <-> LTSA: deployment + disk number     (NFC_A_02 + disk01),
%                  plus the LTSAs of the disks just before and after, so a
%                  session that runs over a disk boundary can still be read.
%
% All four folders are searched recursively:
%   p.tpwsDir (*TPWS*.mat), p.ltsaDir (*.ltsa), p.acceptedDir and
%   p.rejectedDir (*.mat).
%
% Output: struct array, one element per disk, sorted by deployment, disk
% number and letter, with fields
%   key, depRaw, dep, diskNum, letter, label (e.g. NFC_A_02 disk01a)
%   tpws      - full path of the TPWS file
%   tpwsName  - TPWS file name
%   idAcc     - accepted ID file ('' if none)
%   idRej     - rejected ID file ('' if none)
%   idName    - ID file name used for checked\ and checked\inProgress\
%   ltsa      - cell of LTSA file paths (this disk, and its neighbors)
%   ltsaSess  - LTSA session file in p.outDir

for f = {'tpwsDir','acceptedDir','rejectedDir','outDir'}
    if ~isfield(p,f{1}) || isempty(p.(f{1}))
        error('detEdit:acceptReject','Accept/reject review needs p.%s to be set.',f{1});
    end
end

acc = find_files(p.acceptedDir,'*.mat','accepted ID');
rej = find_files(p.rejectedDir,'*.mat','rejected ID');
tpws = find_files(p.tpwsDir,'*TPWS*.mat','TPWS');
if isfield(p,'ltsaDir') && ~isempty(p.ltsaDir)
    ltsa = find_files(p.ltsaDir,'*.ltsa','LTSA');
else
    ltsa = struct('path',{},'name',{},'info',{});
end

% one ID file per key (warn about duplicates)
accMap = key_map(acc,'accepted');
rejMap = key_map(rej,'rejected');
allKeys = unique([keys(accMap), keys(rejMap)]);

tpwsKeys = arrayfun(@(x) x.info.key,tpws,'UniformOutput',false);
ltsaDep = arrayfun(@(x) x.info.dep,ltsa,'UniformOutput',false);
ltsaNum = arrayfun(@(x) x.info.diskNum,ltsa);

files = struct('key',{},'depRaw',{},'dep',{},'diskNum',{},'letter',{},'label',{},...
    'tpws',{},'tpwsName',{},'idAcc',{},'idRej',{},'idName',{},'ltsa',{},'ltsaSess',{});
for iK = 1:length(allKeys)
    key = allKeys{iK};
    idAcc = ''; idRej = '';
    if isKey(accMap,key), idAcc = accMap(key).path; idInfo = accMap(key); end
    if isKey(rejMap,key), idRej = rejMap(key).path; if isempty(idAcc), idInfo = rejMap(key); end, end

    % TPWS: same deployment + disk + letter. If several match, take the one
    % whose name is closest to the ID file name (longest shared start).
    hit = find(strcmp(tpwsKeys,key));
    if isempty(hit)
        fprintf('[WARN] No TPWS file for %s (%s) - skipped.\n',idInfo.name,disk_label(idInfo.info));
        continue
    end
    if length(hit) > 1
        shared = arrayfun(@(h) shared_start(lower(tpws(h).name),lower(idInfo.name)),hit);
        [~,best] = max(shared);
        fprintf('[WARN] %d TPWS files match %s - using %s\n',length(hit),...
            disk_label(idInfo.info),tpws(hit(best)).name);
        hit = hit(best);
    end

    e = struct();
    e.key = key;
    e.depRaw = tpws(hit).info.depRaw;
    e.dep = tpws(hit).info.dep;
    e.diskNum = tpws(hit).info.diskNum;
    e.letter = tpws(hit).info.letter;
    e.label = disk_label(tpws(hit).info);
    e.tpws = tpws(hit).path;
    e.tpwsName = tpws(hit).name;
    e.idAcc = idAcc;
    e.idRej = idRej;
    e.idName = idInfo.name;
    % LTSAs: this disk and the ones before/after it, in name order
    lt = find(strcmp(ltsaDep,e.dep) & abs(ltsaNum - e.diskNum) <= 1);
    [~,o] = sort(lower({ltsa(lt).name}));
    e.ltsa = {ltsa(lt(o)).path};
    e.ltsaSess = fullfile(p.outDir,strrep(e.tpwsName,'TPWS','LTSA'));
    files(end+1) = e; %#ok<AGROW>
end

% sort: deployment, disk number, letter
if ~isempty(files)
    sortKey = arrayfun(@(x) sprintf('%s|%06d|%s',x.dep,x.diskNum,x.letter),files,...
        'UniformOutput',false);
    [~,o] = sort(sortKey);
    files = files(o);
end
fprintf('Accept/reject review: %d disk(s) with accepted/rejected labels and a TPWS file\n',...
    length(files));
end

% ----------
function out = find_files(folder,pattern,what)
% recursive file search; keeps only names with a diskNN tag
out = struct('path',{},'name',{},'info',{});
if ~exist(folder,'dir')
    error('detEdit:acceptReject','%s folder not found: %s',what,folder);
end
[paths,names] = dir_recursive(folder,pattern);
nNoDisk = 0;
for i = 1:length(names)
    info = ar_parse_name(names{i});
    if ~info.ok
        nNoDisk = nNoDisk + 1;
        continue
    end
    out(end+1) = struct('path',paths{i},'name',names{i},'info',info); %#ok<AGROW>
end
if nNoDisk > 0
    fprintf('[WARN] %d %s file(s) in %s have no diskNN in the name - ignored.\n',...
        nNoDisk,what,folder);
end
end

% ----------
function m = key_map(list,what)
m = containers.Map('KeyType','char','ValueType','any');
for i = 1:length(list)
    k = list(i).info.key;
    if isKey(m,k)
        fprintf('[WARN] Two %s files for %s - keeping %s, ignoring %s\n',what,...
            disk_label(list(i).info),m(k).name,list(i).name);
        continue
    end
    m(k) = list(i);
end
end

% ----------
function s = disk_label(info)
s = sprintf('%s disk%02d%s',info.depRaw,info.diskNum,info.letter);
end

% ----------
function n = shared_start(a,b)
L = min(length(a),length(b));
n = find(a(1:L) ~= b(1:L),1) - 1;
if isempty(n), n = L; end
end

% ----------
function [paths,names] = dir_recursive(folder,pattern)
% files matching pattern in folder and all its subfolders
d = dir(fullfile(folder,pattern));
d = d(~[d.isdir]);
names = {d.name};
paths = cellfun(@(n) fullfile(folder,n),names,'UniformOutput',false);
sub = dir(folder);
sub = sub([sub.isdir] & ~ismember({sub.name},{'.','..'}));
for i = 1:length(sub)
    [p2,n2] = dir_recursive(fullfile(folder,sub(i).name),pattern);
    paths = [paths, p2]; %#ok<AGROW>
    names = [names, n2]; %#ok<AGROW>
end
end
