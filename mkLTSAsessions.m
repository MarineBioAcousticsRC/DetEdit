function mkLTSAsessions(userFunc)
% mkLTSAsessions.m
% 2/21/15 version 1.1
% use individual click detections to define session/bout
% get and save ltsa pixel data for each session
% 140310 smw
% clear all

    % MTT = click detection time, MPP = click RL level (dBpp counts)
    % MSN = Nclicks x length snip, MSP = Nclicks x spectra
tic % start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:

% get user input and set up function name
typeInput = exist('userFunc','var');
if typeInput ~= 1
    [userfile,userpath] = uigetfile('*.m',...
        'Select Script with your Data Parameter Settings');
    addpath(userpath) % it adds user folder path to the beggining of the set path
    userFunc = str2func(['@',userfile(1:end-2)]);
end

p = getParams(userFunc,'analysis','mkLTSA');

%% Define subfolder that fit specified iteration
if p.iterationNum > 1
    for id = 2: str2num(p.iterationNum) % iterate id times according to p.iterationNum
        subfolder = ['TPWS',num2str(id)];
        p.tpwsDir = (fullfile(p.tpwsDir,subfolder));
    end
end

%% Check if TPWS file exists (does not look in subdirectories)
% Concatenate parts of file name
if isempty(p.speName)
    detfn = [p.filePrefix,'.*','TPWS',p.iterationNum,'.mat'];
else
    detfn = [p.filePrefix,'.*',p.speName,'.*TPWS',p.iterationNum,'.mat'];
end
% Get a list of all the files in the start directory
fileList = cellstr(ls(p.tpwsDir));
% Find the file name that matches the p.filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Handle Transfer Function
% add in transfer function if desired
if (p.tfSelect > 0)
    [tf,~,~] = getTransfunc(p.filePrefix, p.tfName,p);
else
    tf = 0;
    disp('No TF Applied ');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make LTSA.mat file for each TPWS file found
for iD = 1:length(fileMatchIdx)
    
    %% Find files and load data
    % Find if TPWS file exist
    matchingFile = fileList{fileMatchIdx(iD)};
    detfn = dir(fullfile(p.tpwsDir,matchingFile));
    
    fNameList.TPWS = fullfile(p.tpwsDir,detfn.name);
    A1 = exist(fNameList.TPWS,'file');
    if A1 ~= 2
        disp(['Error: File Does Not Exist: ',fNameList.TPWS])
        return
    end
    
    %%% special case for GOM files (rewording file names)
    % ltsa file for GOM are named with GofMX and without underscore, and tf
    % as well.
    % E.g. GOM_DT_09 -> ltsa is GofMX_DT09
    if strfind(detfn.name,'GOM')
        unscores = strfind(p.filePrefix,'_');
        p.filePrefix(unscores(2)) = '';
        p.filePrefix = regexprep(p.filePrefix,'GOM','GofMX');
    end
    
    % Load detections
    load(fNameList.TPWS)
    
    % test for duplicates in MTT
    ia = []; ic = [];
    [uMTT,ia,ic] = unique(MTT);
    if (length(uMTT) ~= length(MTT))
        disp([' TimeLevel Data NOT UNIQUE - removed:   ', ...
            num2str(length(ic) - length(ia))]);
    end
    [r,c] = size(MTT); %get shape of array
    if (r > c)
        clickTimes = MTT(ia);
        clickLevels = MPP(ia);
    else
        clickTimes = MTT(ia)';
        clickLevels = MPP(ia)';
    end
    
    % Apply tf (if defined) and remove low amplitude detections
    clickLevels = clickLevels + tf;
    ib = find(clickLevels >= p.threshRL);
    disp([' Removed too low:',num2str(length(ia)-length(ib))]);
    clickTimes = clickTimes(ib);
    clickLevels = clickLevels(ib);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Get ltsa file names for a specific site name and deployment number
    d = dir(fullfile(p.ltsaDir,[p.filePrefix,'*.ltsa']));
    fnames = char(d.name);
    nltsas = length(d);
    
    % Load up rawfile start times
    doff = datenum([2000 0 0 0 0 0]);   % convert ltsa time to millenium time
    sTime = []; eTime = []; rfTime = []; hdrStore = [];
    
    if isempty(sTime)
        if nltsas > 0
            sTime = zeros(nltsas,1); eTime = zeros(nltsas,1);
            disp('reading ltsa headers, please be patient ...')
            for k = 1:nltsas
                hdrStore{k} = ioReadLTSAHeader(fullfile(p.ltsaDir,fnames(k,:)));
                sTime(k) = hdrStore{k}.ltsa.start.dnum + doff;  % start time of ltsa files
                eTime(k) = hdrStore{k}.ltsa.end.dnum + doff;    % end time of ltsa files
                rfTime{k} = hdrStore{k}.ltsa.dnumStart + doff; % all rawfiles times for all ltsas
            end
            disp('done reading ltsa headers')
        else
            disp(['No LTSAs found to match wildcard: ', fullfile(p.ltsaDir,[p.filePrefix,'*'])])
            return
        end
    end
    
    % Define LTSA session output file
    lsfn = strrep(detfn.name,'TPWS','LTSA');
    fn2 = fullfile(p.tpwsDir,lsfn);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate bout starts and ends
    
    [nb,eb,sb,bd] = calculate_bouts(clickTimes,p);
    
    % Define bouts if minimum bout duration is given
    if ~isempty(p.minDur)
        minbcount = 1;  %counter for while loop
        while minbcount <= nb
            if bd(minbcount) <= (p.minDur/(60*24))
                sb(minbcount) = sb(minbcount) - ((p.minDur./2)./(60*24)); %subtracts half p.minDur from end time
                eb(minbcount) = eb(minbcount) + ((p.minDur./2)./(60*24)); %adds half p.minDur to end time
            end
            minbcount = minbcount + 1;
        end
        % only if start or end time pass the ltsa time after adding minimum
        % duration
        sb(sb<sTime(1)) = sTime(1);
        eb(eb>eTime(end)) = eTime(end);
    end
    
    %%% Main Loop
    pt = [];
    pwr = [];
    % Loop over the number of bouts (sessions)
    k = 1;
    while (k <= nb)
        % find which ltsa to use and get pwr and pt vectors
        K = [];
        K = find(sTime <= sb(k) & eTime >= eb(k));
        if length(K) > 1
            K = find(sTime(k) <= sb(k) & eTime(k) >= eb(k));
        end
        % find which rawfiles to plot ltsa
        if ~isempty(K) && length(K) == 1
            L = [];
            if eb(k) -sb(k) < p.rawFileDur/ (60*60*24)
                L = find(rfTime{K} >= sb(k),1,'first');
            else
                L = find(rfTime{K} >= sb(k) & rfTime{K} <= eb(k));
            end
            if ~isempty(L)
                if L ~= 1
                    L = [L(1)-1,L]; % get rawfile from before sb(k)
                end
                [pwr{k},pt{k}] = ioGetPwrPtLTSA(p,fnames(K,:),L,K,rfTime,k,hdrStore{K});
            else
                rfT = rfTime{K};
                disp('L is empty')
                pt{k} = [];
                pwr{k} = [];
                disp(['bout start time is ',datestr(sb(k))])
                disp(['bout end time is ',datestr(eb(k))])
            end
        % use the end of one and the beginning of next ltsa   
        elseif isempty(K)   
            disp('K is empty - session spans two LTSAs')
            disp(['bout start time is ',datestr(sb(k))])
            disp(['bout end time is ',datestr(eb(k))])
            
            Ks = [];
            Ks = find(sTime <= sb(k) & eTime >= sb(k));
            Ke = [];
            Ke = find(sTime <= eb(k) & eTime >= eb(k));
            if isempty(Ks) || isempty(Ke)
                disp('Error: Ks or Ke are empty')
                k = k+1;
                continue
            end
            
            Ls = [];
            Ls = find(rfTime{Ks} >= sb(k));
            Le = [];
            Le = find(rfTime{Ke} <= eb(k));
            if ~isempty(Ls)
                if Ls ~= 1
                    Ls = [Ls(1)-1,Ls]; % get rawfile from before sb(k)
                end
                [pwrLs,ptLs] = ioGetPwrPtLTSA(p,fnames(Ks,:),Ls,Ks,rfTime,[],hdrStore{Ks});
            end
            if ~isempty(Le)
                [pwrLe,ptLe] = ioGetPwrPtLTSA(p,fnames(Ke,:),Le,Ke,rfTime,[],hdrStore{Ke});
            end
            
            if isempty(Ls) || isempty(Le)
                disp('Error: Ls or Le are empty ')
                k = k + 1;
                pwr{k} = [];
                pt{k} = [];
                continue
            % combine from end of ltsa with begin of next
            else            
                pwr{k} = [pwrLs pwrLe];
                pt{k} = [ptLs ptLe];
            end
        else
            disp(['K = ',num2str(K')])
            disp(['bout start time is ',datestr(sb(k))])
            disp(['bout end time is ',datestr(eb(k))])
        end
        disp(['Session: ',num2str(k),'  Start: ',datestr(sb(k)),'  End:',datestr(eb(k)),...
            '   Duration: ',num2str(24*60*60*bd(k)),' sec'])
        k = k+1;
    end

    if nb>= 1
        save(fn2,'pwr','pt','-v7.3')
        disp(['Done with file ',fNameList.TPWS])
        tc = toc;
        disp(['Elasped Time : ',num2str(tc),' s'])
    end
    
end

end
