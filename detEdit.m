% detEdit.m
% 2/21/2016 modified for version 1.1
% For Kogia JAH 5/22/15
% Estimate the number of False Detections
% JAH 10-19-2014
% spec2 uses the LTSA for the click spectra
% spec3 uses the TPWS file for the click spectra  JAH 9-26-14
% spcc4 used the TPWS2 file JAH 10-12-14
% 7-7-14 use Simone bouts and Sean Detector JAH
% includes brushing FD, MD in and out of files
% Brushing only works in MATLAB ver 2013b, not 2013a or 2012b
% modified for BW 140308 jah 140320 jah for small ici
% 140311 smw detection editor based on evalSessions.m
clearvars
% utSetDesktopTitle('detEdit');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load default settings
spParamsUser = [];

% Load user input. Has to happen first so you know species.
detEdit_settings


%% Get user input
% and set up file names, if not provided by settings
if ~exist('stn','var')
    stn = input('Enter Project Name (MC GC DT SOCAL): ','s'); % site name
elseif~exist('dpn','var')
    dpn = input('Enter Deployment number + site (01 02 ...): ','s'); % deployment number
elseif~exist('itnum','var')
    itnum = input('Enter Iteration number (1 2 ...): ','s');
elseif~exist('srate','var')
    srate = input('Enter sample rate in kHz (eg 200 or 320): ');
elseif ~exist('sp','var')
    sp = input('Enter Species: Zc Me BWG Md Ko De Po ','s');
elseif ~exist('c4fd','var')
    c4fd = input('Enter Interval to test for False Detections: ') ; %check 4 fd
elseif ~exist('sdir','var')% Get Directory with Detections
    disp('Select Directory with Detections');
    sdir = uigetdir('I:\','Select Directory with Detections');
end

% Get parameter settings worked out between user preferemces, defaults, and
% species-specific settings:
p = sp_setting_defaults(sp,spParamsUser);

gt = gth*60*60;    % gap time in sec
sdn = [stn,dpn];    % site name and deployment number

%% Check if TPWS file exists

detpn = [sdir,'\'];
detfn = [stn,dpn,p.speName,'_TPWS',itnum,'.mat'];
fn = fullfile(detpn,detfn);
if exist(fn,'file') == 0
    fprintf('ERROR: No file named %s exists\n',fn)
    break
end
%% Handle Transfer Function
% add in transfer function if desired
if p.tfSelect > 0
    if ~exist('tfName','var')% user interface to get TF file
        disp('Load Transfer Function File');
        [fname,pname] = uigetfile('I:\Harp_TF\*.tf','Load TF File');
    else % or get it from filename provided in settings
        [pname,fname] = fileparts(tfName);
        fname = [fname,'.tf'];
    end
    tffn = fullfile(pname,fname);
    if strcmp(num2str(fname),'0')
        disp('Cancelled TF File');
        return
    else %give feedback
        disp(['TF File: ',tffn]);
    end
    fid = fopen(tffn);
    [A,count] = fscanf(fid,'%f %f',[2,inf]);
    tffreq = A(1,:);
    tfuppc = A(2,:);
    fclose(fid);
    
    tf = interp1(tffreq,tfuppc,p.tfSelect,'linear','extrap');
    disp(['TF @',num2str(p.tfSelect),' Hz =',num2str(tf)]);
else
    tf = 0;
    disp('No TF Applied');
end

%% Generate FD, TD, and ID files if needed
% detpn = [sdir,['\',stn,'_',spe],'\'];


% false Detection file1
f1pn = detpn;
f1fn = [stn,dpn,p.speName,'_FD',itnum,'.mat'];
fn2 = fullfile(f1pn,f1fn);
A2 = exist(fn2,'file');
if (A2 ~= 2)
    zFD(1,1) = 1;
    save(fn2,'zFD');    % create new FD
    disp('Make new FD file');
end

% Test false Detection file
lf1pn = detpn;
tfn = [stn,dpn,p.speName,'_TD',itnum,'.mat'];
fn6 = fullfile(f1pn,tfn);

% ID File
fnID1 = [stn,dpn,p.speName,'_ID',itnum,'.mat'];
fnID = fullfile(f1pn,fnID1);
A2 = exist(fnID,'file');
if (A2 ~= 2)
    zID = [];
    save(fnID,'zID');    % create new FD
    disp('Make new ID file');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load detections and false detections
% MTT = time MPP = peak-peak % MSN = waveform %MSP = spectra
load(fn,'MTT','MPP')


% if you have more than "maxDetLoad" detections, don't load all spectra and
% time series into memory. You can access them from the disk instead.
% Note: Can't remove duplicates in that case, because matlab won't let you
% select non-contiguous sets from files stored on disk.
nDets = length(MTT);

loadMSP = nDets <= maxDetLoad; % true/false, if you have more detections than
% the maximum load this becomes false.
ic1 = [];
if loadMSP 
    % remove duplicates from MTT (can't do this if too many detections to load into memory).
    [uMTT,ia1,ic1] = unique(MTT);
    if (length(uMTT) ~= length(MTT))
        disp([' TimeLevel Data NOT UNIQUE - removed:   ', ...
            num2str(length(ic1) - length(ia1))]);
    end
    load(fn,'MSP','MSN')
else
    ia1 = [1:length(MTT)]';
end

[r,c] = size(MTT); %get shape of array
if (r > c)
    ct = MTT(ia1);
    cl = MPP(ia1);
else
    ct = MTT(ia1)';
    cl = MPP(ia1)';
end

if specploton && loadMSP
    % if specploton and there aren't too many detections, load spectra 
    csn = MSN(ia1,:); 
    csp = MSP(ia1,:);
else
    disp('No Waveform or Spectra');
end

%% apply tf and remove low amplitude detections
cl = cl + tf;
ib1 = find(cl >= p.threshRL);

if (size(ib1,1) ~= size(cl,1)) && ~loadMSP % catch for case where enforcing
    % min RL threshold on large dataset creates non-continuous indices.
    error('detEdit:RL',['Error: Re-run makeTPWS to enforce your minimum peak to peak RL threshold.\n',...
        'You cannot do it here because you have too many detections to load into memory.\n',...
        sprintf('TPWS minimum RL = %d \ndetEdit minimum RL = %d',min(cl),p.threshRL)])
end

% prune by RL only if spectra & waveforms have been loaded

if specploton && loadMSP
    disp([' Removed too low:',num2str(length(ia1)-length(ib1))]);
    ct = ct(ib1);
    cl = cl(ib1);
    keepers = ia1(ib1);
    
    csn = csn(ib1,:);
    csp = csp(ib1,:);
else
    keepers = ia1;
end

%% Make FD file intersect with MTT
load(fn2)  % false detection times zFD
jFD = []; ia = []; ic = [];
if (~isempty(zFD))
    [jFD,ia,ic] = intersect(MTT,zFD);
    rFD = length(zFD) - length(jFD);
    disp([' Removed ',num2str(rFD),' False Det']);
    if size(jFD,1)<size(jFD,2)
        jFD = jFD';
    end
    zFD = jFD;
    if (isempty(zFD))
        zFD(1,1) = 1;  % keep from blowing up later
    end
    save(fn2,'zFD');
end

%% Make ID file intersect with MTT
load(fnID)  % identified detection times zID
if (~isempty(zID))
    [jID,ia,ic] = intersect(MTT,zID(:,1));
    rID = length(zID) - length(jID);
    disp([' Removed ',num2str(rID),' ID Dets']);
    if size(jID,1)<size(jID,2)
        jID = jID';
    end
    zID = zID(ic,:);
    %     if isempty(zID)
    %         zID = [];  % keep from blowing up later
    %     end
    save(fnID,'zID');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate bout starts and ends 
% TODO: this is calculated in mkLTSA, we should save it there instead of 
% recalculating here but would create backward compatibility issues
% TODO: Should be a function
bFlag = 0;
% find edges (start and end times) of bouts or sessions
dt = diff(ct)*24*60*60; % calculate time between detections in seconds

I = find(dt>gt);  % find start of gaps
sb = [ct(1);ct(I+1)];   % start time of bout
eb = [ct(I);ct(end)];   % end time of bout
dd = ct(end)-ct(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb);      % duration of bout in days

% find bouts > 10 sec long
% bd10 = find(bd > 1 / (60*60*24)); % for Kogia 10 sec
% disp(['Number Bouts : ',num2str(length(bd))])
% limit the length of a bout
blim = p.ltsaMax/24;       % 6 hr bout length limit in days
ib = 1;
while ib <= nb
    bd = (eb - sb);   %duration bout in days
    if (bd(ib) > blim)      % find long bouts
        nadd = ceil(bd(ib)/blim) - 1; % number of bouts to add
        for imove = nb : -1: (ib +1)
            sb(imove+nadd)= sb(imove);
        end
        for iadd = 1 : 1: nadd
            sb(ib+iadd) = sb(ib) + blim*iadd;
        end
        for imove = nb : -1 : ib
            eb(imove+nadd) = eb(imove);
        end
        for iadd = 0 : 1 : (nadd - 1)
            eb(ib+iadd) = sb(ib) + blim*(iadd+1);
        end
        nb = nb + nadd;
        ib = ib + nadd;
    end
    ib = ib + 1;
end

disp(['Number Bouts : ',num2str(nb)])

%% Make LTSA session file
lspn = detpn;
lsfn = [sdn,p.speName,'_LTSA',itnum,'.mat'];
fn5 = fullfile(lspn,lsfn);
A5 = exist(fn5,'file');
if A5 ~= 2
    disp(['Error: LTSA Sessions File Does Not Exist: ',fn5])
    return
else
    disp('Loading LTSA Sessions, please wait ...')
    load(fn5)   % LTSA sessions: pwr and pt structures
    disp('Done Loading LTSA Sessions')
    sltsa = size(pt);
    if (sltsa(2) ~= nb)
        disp('Error: Wrong # LTSA session, REMAKE LTSA ')
        return
    end
end

%% Set up Tests for False Detections
ixfd = (1: c4fd : length(ct));  % selected to test for False Det
A6 = exist(fn6,'file');
if (A6 ~= 2)
    zTD = -1.*ones(nb,4);
    save(fn6,'zTD');    % create new FD
    disp(' Make new TD2 file');
else
    load(fn6)
    if (length(zTD(:,1)) ~= nb)
        disp([' Problem with TD file:',fn6]);
        return
    end
end

%% Compute Spectra Plot Parameters
fimin = 0;% TODO: make this configurable
fimax = srate/2 ; % in kHz 100 or 160 kHz

% set ltsa step size
iPwr = 1;
while isempty(pwr{1,iPwr}) && iPwr<length(pwr)
    iPwr = iPwr+1;
end

if exist(p.dfManual,'var') % allow non-standard ltsa step sizes
    df = p.dfManual;
else
    df = 1000*fimax/(size(pwr{1,iPwr},1)-1);
end

% for LTSA PLOT
f = 1000*fimin:df:1000*fimax;
if p.tfSelect > 0 % tfParams isn't being used...
    tfLTSA = interp1(tffreq,tfuppc,f,'linear','extrap')'; % add to LTSA vector
else
    tfLTSA = zeros(size(f))';
end

if specploton
    % check length of MSP
    inFileMat = matfile(fn);
    if ~loadMSP
        MSP = inFileMat.MSP(1,:);
    end
    smsp2 = size(MSP,2);% 2nd element is num fft points
    ift = 1:smsp2;
    % make frequency vector that matches spectral bins
    fmsp = inFileMat.f;
    if isempty(fmsp)
        fmsp = ((srate/2)/(smsp2-1))*ift - (srate/2)/(smsp2-1);
        fprintf('No freq vector in TPWS file. Using approximation based on sample rate.\n')
    end
    % find the indices that are in the range of interest
    fi = find(fmsp > fimin & fmsp <= fimax);
    fimint = fi(1); fimaxt = fi(end);
    % find index of first bin above min freq boundary
    flowt = find(fmsp > p.fLow,1,'first');
    ft = fmsp(fi);
    fmsp = ft;
    
    % for the PP vs RMS plot
    if (p.tfSelect > 0)
        Ptfpp = interp1(tffreq,tfuppc,fmsp*1000,'linear','extrap');
    else
        Ptfpp = zeros(1,smsp2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kstart = input('Starting Session:  ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause on
disp('Press ''b'' key to go backward ')
disp('Press any other key to go forward')
cc = ' ';  % avoids crash when first bout too short
k = kstart;
yell = [];
blag = 0;
%% Main Loop
% loop over the number of bouts (sessions)
while (k <= nb)
    disp([' BEGIN SESSION: ',num2str(k)]);
    % load in FD, MD and TD each session incase these have been modified
    load(fn2); % brings in zFD
    load(fnID); % bring in zID
    %     load(fn3);  % brings in zMD
    load(fn6); % brings in zTD
    if (length(zTD(1,:)) == 2)
        zTD = [zTD,-1.*ones(length(zTD),2)];
        save(fn6,'zTD');
    end
    % Make PP versus RMS plot for all clicks, if all time series are loaded
    figure(51);clf
    h51 = gca;
        
    if specploton && loadMSP
        xmsp0All = csp + repmat(Ptfpp,size(csp,1),1);
        [xmspAll,im] = max(xmsp0All(:,flowt:fimaxt),[],2); % maximum between flow-100kHz
        
        % calculate peak-to-peak amplitude including transfer function
        xmppAll = cl'-tf+ Ptfpp(im + flowt-1); % vectorized version
        plot(h51,xmspAll,xmppAll,'o','MarkerEdgeColor',[.7,.7,.7],'UserData',ct)
        % hold on;  % keep figure(51) plot with hold on
        title(['Based on ',num2str(length(xmppAll)),' clicks']);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find detections and false detections within this bout (session)
    J = []; JFD =[]; Jtrue = []; XFD = []; JID = [];
    J = find(ct >= sb(k) & ct <= eb(k));
    if specploton && loadMSP % have to load consecutive detections if reading from disk
        J = J(1):J(end);
        csnJ = csn(J,:);
        cspJ = csp(J,:);
    elseif specploton % only load the ones you need for this session
        csnJ = inFileMat.MSN(keepers(J),:);
        cspJ = inFileMat.MSP(keepers(J),:);
    end
    
    % get indices  of test clicks in this session
    XFD = find(ct(ixfd) >= sb(k) & ct(ixfd) <= eb(k));
    zTD(k,1) = length(XFD);
    
    % Test for XFD and strcmp('x or z or w') - if no test points skip
    % x = true, z = false, w = window
    if (isempty(XFD) && (strcmp(cc,'x') || ...
            strcmp(cc,'z') || strcmp(cc,'w')));
        disp(' NO Test Detections, so skip')
        k = k + 1;
        continue
    end
    
    if ~isempty(J) % if there are detection in this session
        t = ct(J); % detection times in this session
        disp([' Detection times:',num2str(length(t))]);
        if (~isempty(XFD))
            xt = ct(ixfd(XFD));  %times to test for False Detection
            xPP = cl(ixfd(XFD));   %amplitude for test False Detection
            disp([' Test False Detection times:',num2str(zTD(k,1))]),
        else
            xt = [];
        end
        RL = cl(J);         % received levels in this session
        nd = length(J);     % number of detections in this session
        
        % get false detection times that intersect with detection times
        K2 = []; % holds false indices
        ff2 = 0;
        tfd = [];
        if (~isempty(zFD)) % get times and indices of false detections
            [tfd,K2,~] = intersect(t,zFD(:,1));
            fdRL = RL(K2);
        end
        if ~isempty(K2) % if this session contains false detections
            ff2 = 1; % set false flag to true
            if specploton
                wavFD = mean(csnJ(K2,:),1); % calculate mean false time series
                specFD = cspJ(K2,:); % get set of false spectra
            end
            disp([' False detections:',num2str(length(K2))])
        else
            ff2 = 0;
            disp(' No False Det')
        end
        
        % get ID'd detection times that intersect with detection times
        K3 = []; % holds Id'd indices
        ff3 = 0; % becomes positive if you have ID's detections in this session
        tID = []; % times of ID'd detections
        IDidx = [];
        if ~isempty(zID)
            [tID,K3,IDidx] = intersect(t,zID(:,1));
            rlID = RL(K3);
        end
        if ~isempty(K3)
            ff3 = 1;
            spCodeSet = zID(IDidx,2); % get ID codes for everything in this session
            specIDs = unique(spCodeSet); % get unique ID codes
            if specploton % get mean spectra for each ID'd type
                wavID = [];
                specID = [];
                for iSpID = 1:length(specIDs)
                    thisSet = spCodeSet == specIDs(iSpID);
                    wavID(iSpID,:) = mean(csnJ(K3(thisSet,:),:),1);
                    specID(iSpID,:) = mean(cspJ(K3(thisSet,:),:),1);
                end
            end
            disp([' ID detections:',num2str(length(K3))])
        else
            ff3 = 0;
            disp(' No ID')
        end
        
        % Calculate indices of detections which are neither false nor ID'd
        JFD = J(K2);
        JID = J(K3);
        JFDandID = union(JFD,JID);
        [Jtrue,iJ,~]= setxor(J,JFDandID); % find all true detections
        trueTimes = ct(Jtrue);% vector of true times in this session

        if specploton
            cspJtrue = cspJ(iJ,:); % true spectra in this session
            csnJtrue = csnJ(iJ,:); % true time series in this session
            wtrue = nanmean(csnJtrue,1); % mean of true spectra in this session
            strue = nanmean(cspJtrue,1); % mean of true time series in this session
        end
        disp([' True Detections: ',num2str(length(trueTimes))])
    else
        disp('Error: no detections between bout start and end')
        return
    end
    dt = diff(t)*24*60*60; % inter-detection interval (IDI) and convert from days to seconds
    
    if ff2 % calculate IDI for false and id'd detections
        dtFD = dt(K2(1:end-1));
    end
    if ff3
        dtID = dt(K3(1:end-1));
    end
    
    disp(['END SESSION: ',num2str(k),' Start: ',datestr(sb(k)),...
        ' End:',datestr(eb(k))])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate number of detections per bin
    [KB,binCX,binT,binC] = ndets_per_bin(t,xt,RL,dt,minNdet,nd);
    % filter empty and low number bins
    if isempty(KB) % not sure what this case does?
        disp(['No bins with at least ',num2str(minNdet),' detections'])
        binT = 0;
        binRL = 0;
        binC = 0;
        k = k + 1;  % go to next
        continue
    end
    if (strcmp(cc,'w') && (zTD(k,2) == 0));
        disp(['Session: ',num2str(k),' # Test Detect Bins: ',...
            num2str(length(binCX)),' but NO False']);
        zTD(k,3) = length(binCX);
        zTD(k,4) = 0;
        save(fn6,'zTD');
        k = k + 1;
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number detection per spectral bin in LTSA
    % make a spectra in figure 50
    PT = pt{1,k};   % LTSA session time vector
    pwr1 = pwr{1,k};  % LTSA power vector
    nbinS = length(PT);
    if (nbinS == 0)
        disp('No LTSA for this Session');
        PT(1) = sb(k) ; PT(2) = eb(k); % start end times for plots
        pwr1(1:length(f)) = ones; % make uniform LTSA
    else
        pwr1 = pwr1((1000*fimin/df)+1:round(1000*fimax/df)+1,:)+ repmat(tfLTSA,1,size(pwr1,2)); % limit plot range LTSA
    end
    durS = PT(end) - PT(1);
    
    if specploton
        % allSPEC = norm_spec(cspJ,flowt,fimint,fimaxt);
        figure(50);clf;set(50,'name','Frequency Spectra')
        h50 = gca;
        figure(52);clf;set(52,'name','Waveform')
        h52 = gca;
        if ~isempty(trueTimes)
            % plot average true click spectrum
            trueSpec = norm_spec_simple(cspJtrue,flowt,fimint,fimaxt);
            plot(h50,ft,trueSpec(fimint:fimaxt),'Linewidth',4)
            hold off;
            % average true click waveform
            plot(h52, wtrue + 2*min(wtrue));
            hold off;
        else
            disp(['No true with at least ',num2str(minNdet),' detections'])
        end
        if ff2   % average false click spec
            SPEC2 = norm_spec_simple(specFD, flowt, fimint, fimaxt);
            % plot average false click spectrum
            figure(50); hold on
            plot(ft,SPEC2(fimint:fimaxt),'r','Linewidth',4)
            hold off
            % plot average false click waveform
            figure(52);hold on
            plot(wavFD + min(wavFD),'r');
            hold off
        end
        if ff3  % average id click spec
            specID_norm = [];
            for iSpec = 1:size(specID,1)
                specID_norm(iSpec,:) = norm_spec_simple(specID(iSpec,:), flowt, fimint, fimaxt);
            end
            % plot average ID'd click spectra
            figure(50); hold on
            hID = plot(ft,specID_norm(:,fimint:fimaxt),'Linewidth',4);
            hold off
            
            % plot average ID'd click waveform(s)
            figure(52); hold on
            hID2 = plot((wavID + 5*rand(size(hID))*min(wavID))');
            
            for iC = 1:length(hID) % set colors
                set(hID(iC),'Color',colorTab(specIDs(iC),:))
                set(hID2(iC),'Color',colorTab(specIDs(iC),:))
            end
            hold off
            
        end
        xlabel(h50,'Frequency (kHz)');
        grid(h50,'on')
        xlim(h50, 'manual');
        ylim(h50,[0 1]);
        xlim(h50,[fmsp(1),fmsp(121)])

        xlabel(h52,'Time (1ms @ 200kHz)');
        ylabel(h52,'Amplitude');
        % for all detections in this session, calculate xmpp and xmsp
        xmsp0 = cspJ + ones(length(J),1) *  Ptfpp; % add transfer fun to session's spectra
        [xmsp,im] = max(xmsp0(:,flowt:fimaxt),[],2);  % maximum value
        xmpp = zeros(length(RL),1);
        for imax = 1 : nd % (over number of detections)
            Pmax = Ptfpp(im(imax) + flowt-1);
            xmpp(imax,1) = RL(imax) - tf + Pmax;
        end
        
        % peakFrkHz = fmsp(im)./1000;
        % Plot  PP versus RMS Plot for this session
        figure(51);set(51,'name','RL pp vs. RL rms')
        hold on
        plot(xmsp,xmpp,'.','UserData',t)% true ones in blue
        if ff2 % false in red
            plot(xmsp(K2),xmpp(K2),'r.','UserData',t(K2))
        end
        if ff3 % ID'd in associated color
            for iC2 = 1:length(specIDs) % set colors
                thisIDset = spCodeSet ==specIDs(iC2);
                hPP = plot(xmsp(K3(thisIDset)),xmpp(K3(thisIDset)),'.','UserData',t(K3(thisIDset)));
                set(hPP,'Color',colorTab(specIDs(iC2),:))
            end
        end
       
        
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots stuff now in figure(201)
    warning('off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(201);clf
    set(201,'name','LTSA and time series')
    hA201 = subplot_layout; % Top panel, Figure 201: Received Level
    axes(hA201(1))
    plot(t,RL,'b.','UserData',t)
    hold on
    % plot(t,RL,'b.','UserData',t)
    if ff2 % plot False detections in red
        plot(tfd,fdRL,'r.','UserData',tfd)
        % disp([' false det plotted:',num2str(length(tfd))])
    end
    if ff3 % plot ID'd detections in associated color
        spCodeSet = zID(IDidx,2); % get species codes for everything in this session
        specIDs = unique(spCodeSet); % get unigue species codes
        for iC2 = 1:length(specIDs) % set colors
            thisIDset = spCodeSet ==specIDs(iC2);
            hRLID = plot(tID(thisIDset),rlID(thisIDset),'.','UserData',tID(thisIDset));
            set(hRLID,'Color',colorTab(specIDs(iC2),:))
        end
    end
    hold off
    axis([PT(1) PT(end) p.rlLow p.rlHi])
    datetick('x',15,'keeplimits')
    grid on
    tstr(1) = {fn};
    tstr(2) = {['Session: ',num2str(k),' Start Time ',...
        datestr(sb(k)),' Detect = ',num2str(nd)]};
    title(tstr);
    ylabel('RL [dB re 1\muPa]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(hA201(2))% middle panel LTSA
    c = (p.ltsaContrast/100) .* pwr1 + p.ltsaBright;
    image(PT,f/1000,c)
    axis([PT(1) PT(end) p.ltsaLims])%v2(4)
    set(gca,'yDir','normal')
    ylabel('Frequency (kHz)')
    datetick('keeplimits')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(hA201(3)) % Bottom panel, Figure 201: Inter-Detection Interval
    % make two copies of dt points for brush
    tdt2 = [];
    dt2 = [];
    ldt = length(dt);
    if ldt > 0
        tdt2 = reshape([t(1:ldt),t((1:ldt)+1)]',2*ldt,1);
        dt2 = reshape([dt,dt]',2*ldt,1);
        
        [AX,H1,H2] = plotyy(tdt2,dt2,binT,binC,'plot','semilogy');
        set(H1,'Marker','.','MarkerFaceColor','b','LineStyle','none','UserData',tdt2)
        set(H2,'Marker','o','MarkerFaceColor','c','LineStyle','none',...
            'Markersize',4.5,'UserData',dt2)
        % Note: plotyy is buggy in 2012b, axis handles work only if called
        % using "axes" and avoid calls to "subplot"
        
        % Do setup for 1st axes
        axes(AX(1))
        axis([PT(1) PT(end) 0 p.dtHi])
        datetick('x',15,'keeplimits')
        Ytick = 0:p.dtHi/10:p.dtHi; % make 0.05 Kogia, 0.2 BW
        set(gca,'YTick',Ytick)
        datetick('x',15,'keeplimits')
        grid on
        ylabel('Time between detections [s]')
        
        % Do setup for 2nd axes
        axes(AX(2))
        axis([PT(1), PT(end), 1, 100])
        datetick('x',15,'keeplimits')
        Ytick2 = [.1 1 10 100 1000 10000];
        set(gca,'YTick',Ytick2)
        ylabel('Det/bin')
        xlabel('Time [GMT]')
        title('Inter-Detection Interval (IDI)')
        grid on
        
        %%% plot false and ID
        axes(AX(1))
        hold on
        if ff2
            plot(AX(1),tfd(2:end),dtFD,'.r','UserData',tfd(2:end))
            % no need to double FD since only the blue points are brush captured
        end
        if ff3 % plot ID'd in associated color
            for iC2 = 1:length(specIDs) % set colors
                thisIDset = spCodeSet ==specIDs(iC2);
                hdtID = plot(AX(1),tID(thisIDset(2:end)),dtID(thisIDset(2:end)),'.','UserData',tID(thisIDset));
                set(hdtID,'Color',colorTab(specIDs(iC2),:))
            end
        end
        hold off
    else
        plot(0,0);
    end
    
    % if you have items brushed in yellow, highlight those on each plot
    if specploton && ~isempty(yell) && ~isempty(csnJ)
        figure(201)
        axes(hA201(1))
        hold on
        plot(t(yell),RL(yell),'ko','MarkerSize',6,'UserData',t(yell));
        hold off;
        
        % for diffs, yell can't exceed length dt, which could happen if you
        % grabbed the last point in the vector, so:
        yellDT = yell(yell<length(dt));
        axes(AX(1))
        hold on
        plot(t(yellDT),dt(yellDT),'ko','MarkerSize',6,'UserData',t(yell));
        hold off
        
        figure(51)
        hold on
        plot(xmsp(yell),xmpp(yell),'ko','MarkerSize',10,...
            'LineWidth',2,'UserData',ct(K2))
        hold off
        
        figure(52); % add click to waveform plot in BLACK
        hold on
        plot(mean(csnJ(yell,:),1)','k');
        hold off
     
        figure(50);  % add click to spec plot in BLACK
        hold on
        cspJy = mean(cspJ(yell,:),1);
        tSPEC = norm_spec_simple(cspJy, flowt,fimint,fimaxt);
        plot(ft,tSPEC(:,fimint:fimaxt),'k','Linewidth',4);
        hold off
    end
    % end of plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    yell = [];
    pause  % wait for user input.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get key stroke
    cc = get(gcf,'CurrentCharacter');
    if strcmp(cc,'u') || strcmp(cc,'g') || strcmp(cc,'y') ||  strcmp(cc,'r');
        % detections were flagged by user
        disp(' Update Display') % Stay on same bout
        % get brushed data and figure out what to do based on color:
        [yell,zFD,zID,bFlag] = brush_color(gca, cc, zFD, zID, colorTab, t);
        
    elseif strcmp(cc,'s') % change time diff scale on bottom plot of 201
        p.dtHi = input(' Update IPI scale (sec):  '); % Set IPI scale
        
    elseif strcmp(cc,'d') % change RL scale on top plot of 201
        p.rlLow = input(' Update RL low (dB):  '); % Set RL low
        p.rlHi = input(' Update RL high (dB):  '); % Set RL high
        
    elseif strcmp(cc,'a')% change LTSA parameters
        p.ltsaContrast = input(sprintf('  Current Contrast %d. Update Contrast:  ',p.ltsaContrast)); 
        p.ltsaBright = input(sprintf('  Current Brightness %d. Update Brightness:  ',p.ltsaBright));
        
    elseif strcmp(cc,'b') % Move backward one bout
        if k ~= 1
            k = k-1;
        end
    elseif strcmp(cc,'f') % assign ALL as false
        disp(['Number of False Detections Added = ',num2str(length(trueTimes))])
        if ~isempty(zID)
            %[newFD,~] = setdiff(t,zID(:,1)); % remove from zID
            [~,iCID] = setdiff(zID(:,1),t); % remove from zID 
            zID = zID(iCID,:);
        end
        newFD = t;
        %end
        zFD = [zFD; newFD; trueTimes]; % Add everything to zFD
        
    elseif strcmp(cc,'t') %assign ALL as true
        [zFD,~] = setdiff(zFD(:,1),t);
        disp(['Remaining False Detections = ',num2str(length(zFD))])
        if ~isempty(zID)
            [~,iCID] = setdiff(zID(:,1),t);
            zID = zID(iCID,:);
        end
        
    elseif strcmp(cc,'j')% jump to non-consecutive session
        prompt = 'Jump to Session: ';
        kjump = input(prompt);
        if (kjump > 0 && kjump < nb)
            k = kjump;
        end
           
    elseif (strcmp(cc,'x') || strcmp(cc,'z') ); % test click for random False Detect
        [zFD,zTD] = test_false_dets(XFD,k,zTD,zFD,xt,xPP,ixfd,ft,trueSpec,flowt,...
            fimint,fimaxt,csn,csp,loadMSP,specploton);
        k = k+1;
       
    elseif (strcmp(cc,'w') && (zTD(k,2) > 0));  % test 5 min window
        % Test 5 min window
        zTD = test_false_bins(k,zTD,xt,xPP,binCX);
        k = k+1;
          
    else
        k = k+1;  % move forward one bout
    end
    
    % after edits, remove duplicate labels and save updated vectors
    if ~isempty(zFD)
        zFD = unique(zFD);
    end
    if ~isempty(zID)
        [~,uniqueID] = unique(zID(:,1));
        zID = zID(uniqueID,:);
    end
    save(fn2,'zFD')
    save(fnID,'zID')    
    save(fn6,'zTD');
    
    % don't end if you used paintbrush on last record
    if (k > nb) && bFlag
        k = nb;
        disp(' Last Record')
        
    end
    bFlag = 0;
end
pause off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make zFD unique
uZFD = [];  ia = []; ic = [];
load(fn2);   % load false detections
[uzFD,ia,ic] = unique(zFD);     % make zFD have unique entries
if (length(ia) ~= length(ic))
    disp([' False Detect NOT UNIQUE - removed:   ', ...
        num2str(length(ic) - length(ia))]);
end
zFD = uzFD;
save(fn2,'zFD');
tfinal = find(zTD(:,1) > 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(['Number of Starting Detections = ',num2str(length(ct)+2)])
disp(' ')
disp(['Number of True Detections = ',num2str(length(ct)-length(zFD)+2)])
disp(' ')
disp(['Number of False Detections = ',num2str(length(zFD)-1)])
disp(' ')
% disp(['Number of Mis-ID Detections = ',num2str(length(zMD(:,1))-1)])
% disp(' ')
disp(['Number of Test Detections & False Detect = ',num2str(sum(zTD(tfinal,:)))])
disp(' ')
disp(['Done with file ',fn])

commandwindow;