function detEdit(userFunc)

% detEdit.m

% Main script to display interface, it takes the user parameter settings
% and plots data in different panels to annotate the data.

% Input
%
%   userFunc - Script user parameter settings.
%       Optional, user wil be prompt to select a scrip if not provided
%
% Examples:
%
% detEdit(@yourDataSettings)
%
% detEdit
clear global
global dPARAMS p dHANDLES fNameList zID zFD zTD
%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:

% get user input and set up function name
typeInput = exist('userFunc','var');
if typeInput ~= 1
    % detfault point to settings folder
    thisPath = mfilename('fullpath');
    
    [userfile,userpath] = uigetfile(fullfile(fileparts(thisPath),...
         'Settings\*.m'),'Select Script with your Data Parameter Settings');
    addpath(userpath) % it adds user folder path to the beggining of the set path
    userFunc = str2func(['@',userfile(1:end-2)]);
end

p = getParams(userFunc,'analysis','detEdit');
p.sizePoints = 8; % the current points
p.sizeBack = 5; % the background points
p.sizeFPR = 8; % the current points
p.colorPoints = [1 .84 0];

%% Define subfolder that fit specified iteration
if p.iterationNum > 1
    for id = 2: str2num(p.iterationNum) % iterate id times according to p.iterationNum
        subfolder = ['TPWS',num2str(id)];
        p.tpwsDir = (fullfile(p.tpwsDir,subfolder));
    end
end

%% Check if TPWS file exists
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
    error(sprintf('No files matching file prefix ''%s'' found!',detfn))
elseif length(fileMatchIdx)>1
    % if more than one match, throw error
    error(sprintf('Multiple TPWS files match the file prefix ''%s''.\n Make the prefix more specific.',detfn))
end

matchingFile = fileList{fileMatchIdx};

%% Handle Transfer Function
% add in transfer function if desired
if p.tfSelect > 0
    [dPARAMS.tf,tffreq,tfuppc] = getTransfunc(p.filePrefix, p.tfName,p);
else
    dPARAMS.tf = 0;
    disp('No TF Applied');
end
%% Generate FD, TD, and ID files if needed
[zFD,zID,fNameList] = buildLabelFiles(matchingFile, p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load detections and false detections
% MTT = time MPP = peak-peak % MSN = waveform %MSP = spectra
fNameList.TPWS = fullfile(p.tpwsDir,matchingFile);
load(fNameList.TPWS,'MTT','MPP')

if isrow(MTT); MTT = MTT'; end
if isrow(MPP); MPP = MPP'; end

% if you have more than "maxDetLoad" detections, don't load all spectra and
% time series into memory. You can access them from the disk instead.
% Note: Can't remove duplicates in that case, because matlab won't let you
% select non-contiguous sets from files stored on disk.
nDets = length(MTT);

if ~isfield(p,'maxDetLoad')
    p.loadMSP = true;
elseif isempty(p.maxDetLoad)
    p.loadMSP = true;
else
    p.loadMSP = nDets <= p.maxDetLoad; % true/false, if you have more detections than
end
% the maximum load this becomes false.
ic1 = [];
if p.loadMSP
    % remove duplicates from MTT (can't do this if too many detections to load into memory).
    [uMTT,ia1,ic1] = unique(MTT);
    if (length(uMTT) ~= length(MTT))
        disp([' TimeLevel Data NOT UNIQUE - removed:   ', ...
            num2str(length(ic1) - length(ia1))]);
    end
    load(fNameList.TPWS,'MSP','MSN')
else
    ia1 = [1:length(MTT)]';
end

[r,c] = size(MTT); %get shape of array
if (r > c)
    dPARAMS.clickTimes = MTT(ia1);
    dPARAMS.clickLevels = MPP(ia1);
else
    dPARAMS.clickTimes = MTT(ia1)';
    dPARAMS.clickLevels = MPP(ia1)';
end

if p.specploton && p.loadMSP
    % if p.specploton and there aren't too many detections, load spectra
    dPARAMS.csn = MSN(ia1,:);
    dPARAMS.csp = MSP(ia1,:);
else
    disp('Number of detections exceeds max for loading');
end

%% Apply tf and remove low amplitude detections
dPARAMS.clickLevels = dPARAMS.clickLevels + dPARAMS.tf;
ib1 = find(dPARAMS.clickLevels >= p.threshRL);

if (size(ib1,1) ~= size(dPARAMS.clickLevels,1)) && ~p.loadMSP % catch for case where enforcing
    % min RL threshold on large dataset creates non-continuous indices.
    error('detEdit:RL',['Error: Re-run makeTPWS to enforce your minimum peak to peak RL threshold.\n',...
        'You cannot do it here because you have too many detections to load into memory.\n',...
        sprintf('TPWS minimum RL = %d \ndetEdit minimum RL = %d',min(dPARAMS.clickLevels),p.threshRL)])
end

if (size(ib1,1) == 0)
    % min RL threshold excludes all detections
    error('detEdit:RL',['Error: No detections meet the minimum peak to peak RL threshold.\n',...
        sprintf('TPWS maximum RL = %d \ndetEdit minimum RL = %d',max(dPARAMS.clickLevels),p.threshRL)])
end

% prune by RL only if spectra & waveforms have been loaded
if p.specploton && p.loadMSP
    disp([' Removed too low:',num2str(length(ia1)-length(ib1))]);
    dPARAMS.clickTimes = dPARAMS.clickTimes(ib1);
    dPARAMS.clickLevels = dPARAMS.clickLevels(ib1);
    dPARAMS.keepers = ia1(ib1);
    
    dPARAMS.csn = dPARAMS.csn(ib1,:);
    dPARAMS.csp = dPARAMS.csp(ib1,:);
else
    dPARAMS.keepers = ia1;
end

%% Make FD file intersect with MTT
load(fNameList.FD)  % false detection times zFD
zFD = rmUnmatchedDetections(MTT, zFD);
save(fNameList.FD,'zFD');

% Make ID file intersect with MTT
load(fNameList.ID)  % identified detection times zID
zID = rmUnmatchedDetections(MTT, zID);
% if exist('labels','var')
%     % a label struct was fount in ID file. Add it to p. This will overwrite
%     % whatever was in mySpID before.
%     p.mySpID = labels;
% elseif exist('mySpID','var')
%     p.mySpID = mySpID;
% end
save(fNameList.ID,'zID','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate bout starts and ends
% TODO: this is calculated in mkLTSA, we should save it there instead of
% recalculating here but would create backward compatibility issues

[dPARAMS.nb,dPARAMS.eb,dPARAMS.sb,dPARAMS.bd] = calculate_bouts(dPARAMS.clickTimes,p);
dPARAMS.bFlag = 0;

%% Make LTSA session file
lsfn = strrep(matchingFile,'TPWS','LTSA');
fNameList.LTSA = fullfile(p.tpwsDir,lsfn);
Altsa = exist(fNameList.LTSA,'file');
if Altsa ~= 2
    disp(['Error: LTSA Sessions File Does Not Exist: ',fNameList.LTSA])
    return
else
    disp('Loading LTSA Sessions, please wait ...')
    load(fNameList.LTSA)   % LTSA sessions: pwr and pt structures
    dPARAMS.pwr = pwr;
    dPARAMS.pt = pt;
    sltsa = size(dPARAMS.pt);
    
    if (sltsa(2) ~= dPARAMS.nb)
        disp(['Error: Number of LTSA sessions calculated here doesn''t match ',...
            'input LTSA file. Check ltsaMax parameter.'])
        dPARAMS.pwr = pwr(1:dPARAMS.nb);
        dPARAMS.pt = pt(1:dPARAMS.nb);
        % return
    end
    
    clear pt pwr
    
    disp('Done Loading LTSA Sessions')


end

%% Set up Tests for False Detections
% % The false positive estimate tool picks every Nth click to test. If you
% % have false positives in zFD, you can pick out only the true ones to
% % determine which click indices to look at (that happens just below),
% % but if the user then adds or removes anything from zFD, the indices won't
% % adjust. Provide a warning to tell the user there might be an issue.
if ~isempty(zFD)
    disp(strcat('WARNING: This dataset contains false-flagged detections.  ', ...
        'Remove them using modDet prior to estimating false positive rate.'))
end

[~,trueClickIDx] = setdiff(dPARAMS.clickTimes, zFD);
% select clicks for each label to test for False Det
ixfd = {};
dPARAMS.testClickIdx = {};
if ~isempty(zID)
for i = 1:length(p.mySpID)
    thisLabTimes = zID(zID(:,2)==p.mySpID(i).zID_Label,1);
    [C iCT iLabT] = intersect(dPARAMS.clickTimes(trueClickIDx),thisLabTimes);
    nClick = size(iCT,1);
    if nClick >= p.c4fd % test p.c4fd clicks
        ixfd{i} = sort(datasample(iCT,p.c4fd,1,'Replace',false)); 
    elseif nClick < p.c4fd % if there aren't enough to reach p.c4fd, test all clicks
        ixfd{i} = iCT;
    end
    dPARAMS.testClickIdx{i} = trueClickIDx(ixfd{i});
end
end
% 
% A6 = exist(fNameList.TD,'file');
% if (A6 ~= 2)
%     zTD = {};
%     cMat = {};
%     fpfnTD = {};
%     for j = 1:dPARAMS.nb
%         for i = 1:length(p.mySpID)
%             zTD{j,1} = j;
%             zTD{j,i+1} = -1.*ones(1,6);
%             fpfnTD{1,i} = [];
%             for k = 1:length(p.mySpID)+1
%                cMat{i,k} = zeros(1,2); 
%             end
%         end
%     end
% 
%     save(fNameList.TD,'zTD','cMat','fpfnTD');    % create new TD
%     disp(' Make new TD file');
% else
%     load(fNameList.TD)
%     if (~exist('zTD') || ~exist('cMat') || ~exist('fpfnTD')) || (size(zTD,1) ~= dPARAMS.nb) || size(zTD,2) < length(p.mySpID)+1
%         disp([' Problem with existing TD file: ',fNameList.TD]);
%         return
%     end
% end

%% Set up False Positive & False Negative Bin Tests

A6 = exist(fNameList.TD,'file');
savep = p;
if (A6 ~= 2)
    cMat = {};
    fpfnTD = {};
    for j = 1:dPARAMS.nb
        for i = 1:length(p.mySpID)
            zTD{j,1} = j;
            zTD{j,i+1} = -1.*ones(1,6);
            fpfnTD{1,i} = [];
            for k = 1:length(p.mySpID)+1
               cMat{i,k} = zeros(1,2); 
            end
        end
    end

    save(fNameList.TD,'zTD','cMat','fpfnTD');    % create new TD
    disp(' Make new TD file');
else
    load(fNameList.TD)
    if (~exist('zTD') || ~exist('cMat') || ~exist('fpfnTD')) || (size(zTD,1) ~= dPARAMS.nb) || size(zTD,2) < length(p.mySpID)+1
        disp([' Problem with existing TD file: ',fNameList.TD]);
        return
    end
end
p = savep;
global cMat fpfnTD

% divide entire TPWS into bins aligned with cluster_bins 
dPARAMS.binTimes = (floor(MTT(1)):datenum([0,0,0,0,p.binDur,0]):ceil(MTT(end)))';
dPARAMS.ftb = {};

% divy all clicks into bins
% remove clicks below necessary MPP threshold (w/e was set in cluster bins)
dPARAMS.clickTimesforBinning = dPARAMS.clickTimes(MPP>=120);

[Ntot,~] = histcounts(dPARAMS.clickTimesforBinning,dPARAMS.binTimes);
%remove bins without sufficient # clicks
binsWithClicks = find(Ntot>p.minClicks);

if ~isempty(zID)
for i = 1:length(p.mySpID) % for each label
    ftb = [];
    % divy clicks with this label into bins
    thisLabTimes = zID(zID(:,2)==p.mySpID(i).zID_Label,1);
    [N,~] = histcounts(thisLabTimes,dPARAMS.binTimes);
    thisClickBins = find(N>0);
    notThisClickBins = setdiff(binsWithClicks,thisClickBins);
    
    
    % randomly select test bins across TPWS from those bins which DO and DO
    % NOT contain clicks with this label
    if p.nTestBins < size(thisClickBins,2)
        ftb = sort(datasample(thisClickBins,p.nTestBins,2,'Replace',false));
    else
        ftb = thisClickBins;
    end
        if p.nTestBins < size(notThisClickBins,2)
        ftb = [ftb,datasample(notThisClickBins,p.nTestBins,2,'Replace',false)];
    else
        ftb = [ftb,notThisClickBins];
        end
    dPARAMS.ftb{1,i} = sort(ftb);
end
end
%% Set up Label Certainty Evaluation 

% divide entire TPWS into bins aligned with cluster_bins
dPARAMS.binTimes = (floor(MTT(1)):datenum([0,0,0,0,p.binDur,0]):ceil(MTT(end)))';
dPARAMS.ixtb = {};

if ~isempty(zID)
for i = 1:length(p.mySpID) %for each label
    % divvy labeled clicks into bins
    thisLabTimes = zID(zID(:,2)==p.mySpID(i).zID_Label,1);
    [N,~] = histcounts(thisLabTimes,dPARAMS.binTimes);
    goodBins = find(N>5);
    
    % randomly select test bins across TPWS
    if p.nTestBins < size(goodBins,2)
        dPARAMS.ixtb{1,i} = sort(datasample(goodBins,p.nTestBins,2,'Replace',false));
    else
        %fprintf(['WARNING: Not enough bins for label ' num2str(p.mySpID(i).zID_Label),' to meet user setting for Label Certainty Evaluation, using all available bins\n']);
        dPARAMS.ixtb{1,i} = goodBins;
    end
    
end
end

%% Compute Spectra Plot Parameters
% max and min for LTSA frequency
dPARAMS.fiminLTSA = 0;% TODO: make this configurable
dPARAMS.fimaxLTSA = p.sampleRate/2 ; % in kHz 100 or 160 kHz

% set ltsa step size
iPwr = 1;
while isempty(dPARAMS.pwr{1,iPwr}) && iPwr<length(dPARAMS.pwr)
    iPwr = iPwr+1;
end

% ToDo: Seems like LTSA parameters (step size and frequency bins) could
% be calculated from LTSA directly, especially because there are
% inconsistent step sizes in some LTSAs.
if any(strcmp('dfManual',fieldnames(p)))&& ~isempty(p.dfManual)
    % allow non-standard ltsa step sizes
    dPARAMS.df = p.dfManual;
else
    dPARAMS.df = 1000*dPARAMS.fimaxLTSA/(size(dPARAMS.pwr{1,iPwr},1)-1);
end

% for LTSA PLOT
dPARAMS.f = 1000*dPARAMS.fiminLTSA:dPARAMS.df:1000*dPARAMS.fimaxLTSA;
if p.tfSelect > 0 % tfParams isn't being used...
    dPARAMS.tfLTSA = interp1(tffreq,tfuppc,dPARAMS.f,'linear','extrap')'; % add to LTSA vector
else
    dPARAMS.tfLTSA = zeros(size(dPARAMS.f))';
end

% check length of MSP
dPARAMS.inFileMat = matfile(fNameList.TPWS);
if ~p.loadMSP
    MSP = dPARAMS.inFileMat.MSP(1,:);
end
smsp2 = size(MSP,2);% 2nd element is num fft points
ift = 1:smsp2;
% make frequency vector that matches spectral bins
if any(strcmp('f',fieldnames(dPARAMS.inFileMat)))
    dPARAMS.fmsp = dPARAMS.inFileMat.f;
else
    dPARAMS.fmsp = [];
end
if isempty(dPARAMS.fmsp)
    dPARAMS.fmsp = ((p.sampleRate/2)/(smsp2-1))*ift - (p.sampleRate/2)/(smsp2-1);
    fprintf('No freq vector in TPWS file. Using approximation based on sample rate.\n')
end
% find the indices that are in the range of interest
fi = find(dPARAMS.fmsp > p.fLow &...
    dPARAMS.fmsp <= p.fHi);
dPARAMS.fimint = fi(1); dPARAMS.fimaxt = fi(end);
dPARAMS.ft = dPARAMS.fmsp(fi);

% for the PP vs RMS plot
if (p.tfSelect > 0)
    dPARAMS.Ptfpp = interp1(tffreq,tfuppc,dPARAMS.fmsp*1000,'linear','extrap');
else
    dPARAMS.Ptfpp = zeros(1,smsp2);
end
if ~isrow(dPARAMS.Ptfpp); dPARAMS.Ptfpp = dPARAMS.Ptfpp'; end

compute_transformations
if p.autoFalse
    % apply automatic false thresholds based on frequency and RMS
    % amplitude
    apply_auto_thresh
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPARAMS.k = input('Starting Session:  ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press ''b'' key to go backward ')
disp('Press any other key to go forward')
dPARAMS.cc = ' ';  % avoids crash when first bout too short
dPARAMS.yell = [];
dPARAMS.blag = 0;

% if (size(zTD,2) == 2) % this seems to patch on extra columns
%     % to old zTD matrices that maybe only had the first two. Probably
%     % only needed for backward compatibility
%     zTD = [zTD,-1.*ones(length(zTD),2)];
%     save(fNameList.TD,'zTD');
% end

% initialize ID toggles to on
dPARAMS.NoLabel_Toggle = 'on';
dPARAMS.FD_Toggle = 'on';
[idToggleInit{1:size(p.colorTab,1),1}] = deal('on');
dPARAMS.ID_Toggle = idToggleInit;

% Check if LTSA plot exists, is so, don't reset position
% if ishghandle(201)   
%    existLTSA = 1;
% else
%    existLTSA = 0; 
% end
dHANDLES.LTSAfig = figure(201); colormap(dHANDLES.LTSAfig, jet)
set(dHANDLES.LTSAfig,'name','LTSA and time series',...
    'KeyPressFcn',{@keyAction})
% if ~existLTSA
%     defaultPosLTSA=[.30,0,.55,1];
%     set(dHANDLES.LTSAfig,'Units','normalized')
%     set(dHANDLES.LTSAfig,'Position',defaultPosLTSA)
% end

dHANDLES.hbLTSA = brush(dHANDLES.LTSAfig);
set(dHANDLES.hbLTSA,'Color',[1,1,0],'Enable','off'); % light yellow [.9290 .6940 .1250]
set(dHANDLES.hbLTSA,'ActionPostCallback',{@brushOff,dHANDLES.LTSAfig})
% 
dHANDLES.RMSvPPfig = figure(51);
set(dHANDLES.RMSvPPfig,'name','RL pp vs. RL rms',...
    'KeyPressFcn',@keyAction)
dHANDLES.hbRMSvPP = brush(dHANDLES.RMSvPPfig);
set(dHANDLES.hbRMSvPP,'Color',[1,1,0],'Enable','off'); % light yellow [.9290 .6940 .1250]
set(dHANDLES.hbRMSvPP,'ActionPostCallback',{@brushOff,dHANDLES.RMSvPPfig})

dHANDLES.RMSvFreqfig = figure(53);
set(dHANDLES.RMSvFreqfig,'name','RL rms vs. Peak freq.',...
    'KeyPressFcn',@keyAction)
dHANDLES.hbRMSvFreq = brush(dHANDLES.RMSvFreqfig);
set(dHANDLES.hbRMSvFreq,'Color',[1,1,0],'Enable','off'); % light yellow [.9290 .6940 .1250]
set(dHANDLES.hbRMSvFreq,'ActionPostCallback',{@brushOff,dHANDLES.RMSvFreqfig})

dHANDLES.spectrafig = figure(50);
set(dHANDLES.spectrafig,'name','Frequency Spectra',...
        'KeyPressFcn',@keyAction)

dHANDLES.wavefig = figure(52);
set(dHANDLES.wavefig,'name','Waveform',...
        'KeyPressFcn',@keyAction)


if exist('plotaxes','var')
    dPARAMS.plotaxes = plotaxes;
end


% figures out a max for spectral plot
if p.threshHiFreq ~=0 %any(strcmp('threshHiFreq',fieldnames(p)))
    dPARAMS.ymax = p.threshHiFreq + 1;
else
    dPARAMS.ymax = dPARAMS.fmsp(end);  % yaxis max of plot 53 (Default)
end


%% Main Loop
% loop over the number of bouts (sessions)
boutMotion % run it once to set up first view. After that it will be 
% run after key press only.

% Display ID legend:
% (note since it's not in the loop, if people close it,
% it won't come back in this session.
% put it at end so it can reference an existing figure for default colors
ID_legend
set(dHANDLES.hID,'name','Legend',...
        'KeyPressFcn',@keyAction)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     dPARAMS.onerun = dPARAMS.onerun+1;
%     % get key stroke
%     cc = get(gcf,'CurrentCharacter');
%     
%    
% pause off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
