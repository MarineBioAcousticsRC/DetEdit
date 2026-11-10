function mkLTSAsessions_AR(p)

% mkLTSAsessions_AR
%
% mkLTSAsessions for accept/reject review (p.plotAcceptReject = true),
% called by mkLTSAsessions.
%
% Makes an LTSA session file for every disk that has an accepted and/or
% rejected ID file, across all deployments found under p.ltsaDir:
%   - disks are matched by deployment + disk (ar_build_index), not by name
%   - sessions are built only around clicks with an accepted or rejected
%     label (ar_bout_times + ar_calculate_bouts), the same way detEdit
%     builds them
%   - session files go to p.outDir (TPWS name with TPWS -> LTSA)
%   - disks whose session file already exists in p.outDir are skipped, so
%     you can stop and re-run
%
% Clicks are read and thresholded the same way as in mkLTSAsessions
% (unique MTT, + transfer function, >= p.threshRL), but only MTT and MPP
% are loaded from each TPWS file.

tStart = tic;
if ~isfield(p,'ltsaDir') || isempty(p.ltsaDir)
    error('detEdit:acceptReject','Set ltsaDir (folder with the LTSA files) in your settings.');
end
if ~exist(p.outDir,'dir')
    mkdir(p.outDir);
end
files = ar_build_index(p);

doff = datenum([2000 0 0 0 0 0]);   % convert ltsa time to millenium time
hdrCache = containers.Map('KeyType','char','ValueType','any');
tfCache = containers.Map('KeyType','char','ValueType','double');
nDone = 0; nSkip = 0; nFail = 0;

for iF = 1:length(files)
    e = files(iF);
    fprintf('\n[%d/%d] %s   %s\n',iF,length(files),e.label,e.tpwsName);
    if exist(e.ltsaSess,'file') == 2
        fprintf('  Session file already exists - skipped: %s\n',e.ltsaSess);
        nSkip = nSkip + 1;
        continue
    end
    if isempty(e.ltsa)
        fprintf('  [WARN] No LTSA file for %s disk%02d in %s - skipped.\n',...
            e.depRaw,e.diskNum,p.ltsaDir);
        nFail = nFail + 1;
        continue
    end

    % transfer function (per deployment)
    if p.tfSelect > 0
        if ~isKey(tfCache,e.dep)
            [tfCache(e.dep),~,~] = getTransfunc(e.depRaw,p.tfName,p);
        end
        tf = tfCache(e.dep);
    else
        tf = 0;
    end

    % clicks: same as mkLTSAsessions
    T = load(e.tpws,'MTT','MPP');
    [~,ia] = unique(T.MTT);
    clickTimes = T.MTT(ia); clickTimes = clickTimes(:);
    clickLevels = T.MPP(ia); clickLevels = clickLevels(:) + tf;
    clear T
    boutTimes = ar_bout_times(clickTimes,clickLevels,p.threshRL,e);
    fprintf('  %d labeled clicks above %g dB define the sessions\n',...
        length(boutTimes),p.threshRL);
    if isempty(boutTimes)
        fprintf('  [WARN] No accepted/rejected click above threshRL - no sessions, skipped.\n');
        nFail = nFail + 1;
        continue
    end

    % LTSA headers (cached: neighboring disks share LTSAs)
    nl = length(e.ltsa);
    hdrStore = cell(1,nl); rfTime = cell(1,nl);
    sTime = zeros(nl,1); eTime = zeros(nl,1);
    for k = 1:nl
        if ~isKey(hdrCache,e.ltsa{k})
            fprintf('  Reading LTSA header %s\n',e.ltsa{k});
            hdrCache(e.ltsa{k}) = ioReadLTSAHeader(e.ltsa{k});
        end
        hdrStore{k} = hdrCache(e.ltsa{k});
        sTime(k) = hdrStore{k}.ltsa.start.dnum + doff;
        eTime(k) = hdrStore{k}.ltsa.end.dnum + doff;
        rfTime{k} = hdrStore{k}.ltsa.dnumStart + doff;
    end
    [sTime,o] = sort(sTime);   % time order
    eTime = eTime(o); hdrStore = hdrStore(o); rfTime = rfTime(o);
    fnames = e.ltsa(o);

    [nb,eb,sb,bd] = ar_calculate_bouts(boutTimes,p);
    [pwr,pt] = session_pwr(p,nb,eb,sb,bd,sTime,eTime,rfTime,hdrStore,fnames);

    if nb >= 1
        save(e.ltsaSess,'pwr','pt','-v7.3')
        fprintf('  Saved %d session(s): %s\n',nb,e.ltsaSess);
        nDone = nDone + 1;
    else
        fprintf('  [WARN] No sessions - nothing saved.\n');
        nFail = nFail + 1;
    end
end

fprintf('\nDone: %d saved, %d already there, %d skipped with a warning. %.0f s\n',...
    nDone,nSkip,nFail,toc(tStart));
end

% ----------
function [pwr,pt] = session_pwr(p,nb,eb,sb,bd,sTime,eTime,rfTime,hdrStore,fnames)
% LTSA pixels for each session: same steps as the loop in mkLTSAsessions
pwr = {}; pt = {};

% (p.minDur is not applied here: short sessions are already widened by
% ar_calculate_bouts, the same way detEdit does, so the LTSA lines up)

k = 1;
while (k <= nb)
    % find which ltsa to use and get pwr and pt vectors
    K = find(sTime <= sb(k) & eTime >= eb(k));
    if length(K) > 1
        K = K(1);
    end
    if ~isempty(K)
        if eb(k) - sb(k) < p.rawFileDur/(60*60*24)
            L = find(rfTime{K} >= sb(k),1,'first');
        else
            L = find(rfTime{K} >= sb(k) & rfTime{K} <= eb(k));
        end
        if ~isempty(L)
            if L(1) ~= 1
                L = [L(1)-1,L(:)']; % get rawfile from before sb(k)
            end
            [pwr{k},pt{k}] = ioGetPwrPtLTSA(p,fnames{K},L,K,rfTime,k,hdrStore{K});
        else
            disp('L is empty')
            pt{k} = []; pwr{k} = [];
            disp(['bout start time is ',datestr(sb(k))])
            disp(['bout end time is ',datestr(eb(k))])
        end
    else
        % use the end of one and the beginning of next ltsa
        disp('K is empty - session spans two LTSAs')
        disp(['bout start time is ',datestr(sb(k))])
        disp(['bout end time is ',datestr(eb(k))])
        Ks = find(sTime <= sb(k) & eTime >= sb(k),1);
        Ke = find(sTime <= eb(k) & eTime >= eb(k),1);
        if isempty(Ks) || isempty(Ke)
            disp('Error: Ks or Ke are empty')
            pwr{k} = []; pt{k} = [];
            k = k + 1;
            continue
        end
        Ls = find(rfTime{Ks} >= sb(k));
        Le = find(rfTime{Ke} <= eb(k));
        if isempty(Ls) || isempty(Le)
            disp('Error: Ls or Le are empty ')
            pwr{k} = []; pt{k} = [];
            k = k + 1;
            continue
        end
        if Ls(1) ~= 1
            Ls = [Ls(1)-1,Ls(:)'];
        end
        [pwrLs,ptLs] = ioGetPwrPtLTSA(p,fnames{Ks},Ls,Ks,rfTime,[],hdrStore{Ks});
        [pwrLe,ptLe] = ioGetPwrPtLTSA(p,fnames{Ke},Le,Ke,rfTime,[],hdrStore{Ke});
        pwr{k} = [pwrLs pwrLe];
        pt{k} = [ptLs ptLe];
    end
    disp(['Session: ',num2str(k),'  Start: ',datestr(sb(k)),'  End:',datestr(eb(k)),...
        '   Duration: ',num2str(24*60*60*bd(k)),' sec'])
    k = k + 1;
end
% one cell per session, even if the last ones are empty
if length(pwr) < nb
    pwr{nb} = []; pt{nb} = [];
end
end
