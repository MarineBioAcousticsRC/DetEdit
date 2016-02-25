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
clear all
%utSetDesktopTitle('detEdit');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some parameters
gth = .5;    % gap time in hrs between sessions
gt = gth*60*60;    % gap time in sec
contrast = 282; bright = 50; % for LTSA
rll = 110; rlh = 170; % PP plot window
specploton = 1; %1 = yes spec plot 0 = no spec plot
ltsamax = 6; % length of ltsa window
%df = 100;   LTSA step size Now caLculated from srate and pwr1 (#ltsa sample)
% get user input and set up file names
stn = input('Enter Project Name (MC GC DT SOCAL): ','s'); % site name
dpn = input('Enter Deployment number + site (01 02 ...): ','s'); % deployment number
itnum = input('Enter Iteration number (1 2 ...): ','s');
srate = input('Enter sample rate in kHz (eg 200 or 320): ');
sp = input('Enter Species: Zc Me BWG Md Ko De Po ','s');
if (strcmp(sp,'Ko') || strcmp(sp,'k'))
    specchar = 'K'; %Simone abbreviation for species
    spe = 'Kogia';  tfselect = 80000; % freq used for transfer function
    dl = 0.5;  % scale for IPI display in sec
    flow = 70;   % 70 kHz boundary for spec plot
    thres = 116; % dB threshold
elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
    specchar = 'Z'; %Simone abbreviations for BW species
    spe = 'Cuviers';  tfselect = 40200; % freq used for transfer function
    dl = 1.0; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 121; % dB threshold
    contrast = 200; bright = 30; % for LTSA
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    specchar = 'M'; %Simone abbreviations for BW species
    spe = 'Gervais';  tfselect = 40200; % freq used for transfer function
    dl = 1.0; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 121; % dB threshold
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    specchar = 'G'; %Simone abbreviations for BW species
    spe = 'BWG';  tfselect = 40200; % freq used for transfer function
    dl = 1.0; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 121; % dB threshold
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    specchar = 'D'; %Simone abbreviations for BW species
    spe = 'BW31';  tfselect = 40200; % freq used for transfer function
    dl = 1.0; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 121; % dB threshold
elseif (strcmp(sp,'De') || strcmp(sp,'de'))
    spe = 'Delphin';  tfselect = 0; % already in dB no correction
    dl = 0.5; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 136.9; % dB threshold
    rll = thres - 6.9; rlh = 190;
elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
    spe = 'Porpoise';  tfselect = 0; % already in dB no correction
    dl = 0.5; % scale for IPI display in sec
    flow = 25;   % 25 kHz boundary for spec plot
    thres = 100; % dB threshold
    rll = thres - 5; rlh = 190;
    contrast = 250; bright = 100; % for LTSA
elseif (strcmp(sp,'MFA') || strcmp(sp,'mfa'))
    spe = 'MFA';  tfselect = 4000; % freq used for transfer function
    ltsamax = .5; % length of ltsa window
    thres = 80;
    dl = 2; % scale for IPI display in sec
    flow = 2;   % 2 kHz boundary for spec plot
    rll = thres - 5; rlh = 180;
    ltsamax = 0.5;
    df = 10;    % LTSA step size in 10 [Hz] bins
else
    disp(' Bad Species type')
    return
end
c4fd = input('Enter Interval to test for False Detections: ') ; %check 4 fd
sdn = [stn,dpn];    % site name and deployment number
% user interface to get TF file
disp('Load Transfer Function File');
[fname,pname] = uigetfile('I:\Harp_TF\*.tf','Load TF File');
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
if (tfselect > 0)
    tf = interp1(tffreq,tfuppc,tfselect,'linear','extrap');
    disp(['TF @',num2str(tfselect),' Hz =',num2str(tf)]);
else
    tf = 0;
    disp('No TF Applied');
end
% Get Directory with Detections
disp('Select Directory with Detections');
sdir = uigetdir('I:\','Select Directory with Detections');
%
% detpn = [sdir,['\',stn,'_',spe],'\'];
detpn = [sdir,'\'];
detfn = [stn,dpn,'_',spe,'_TPWS',itnum,'.mat'];
fn = fullfile(detpn,detfn);
% false Detection file1
f1pn = detpn;
f1fn = [stn,dpn,'_',spe,'_FD',itnum,'.mat'];
fn2 = fullfile(f1pn,f1fn);
A2 = exist(fn2,'file');
if (A2 ~= 2)
    zFD(1,1) = 1;
    save(fn2,'zFD');    % create new FD
    disp([' Make new FD file']);
end
% Test false Detection file
lf1pn = detpn;
tfn = [stn,dpn,'_',spe,'_TD',itnum,'.mat'];
fn6 = fullfile(f1pn,tfn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load detections and false detections
% MTT = time MPP= peak-peak % MSN = waveform %MSP = spectra
load(fn)
% test that MTT is unique
ia = []; ic = [];
[uMTT,ia,ic] = unique(MTT);
if (length(uMTT) ~= length(MTT))
    disp([' TimeLevel Data NOT UNIQUE - removed:   ', ...
        num2str(length(ic) - length(ia))]);
end
[r,c] = size(MTT); %get shape of array
if (r > c)
    ct = MTT(ia);
    cl = MPP(ia);
else
    ct = MTT(ia)';
    cl = MPP(ia)';
end
if (specploton == 1) % only if spec plot
    if (strcmp(spe,'Delphin'))% set length of waveform plot to 1 ms
        csn = MSN(ia,100:100+(srate+2));
    else
        csn = MSN(ia,:);
    end
    csp = MSP(ia,:);
else
    disp('No Waveform or Spectra');
end
% apply tf and remove low amplitude
cl = cl + tf;
ib = find(cl >= thres);
disp([' Removed too low:',num2str(length(ia)-length(ib))]);
ct = ct(ib);
cl = cl(ib);
if (specploton == 1) % do this only if spec plot
    csn = csn(ib,:);
    csp = csp(ib,:);
end
% Make FD file intersect with MTT
load(fn2)  % false detection times zFD
jFD = []; ia = []; ic = [];
if (~isempty(zFD))
    [jFD,ia,ic] = intersect(MTT,zFD);
    rFD = length(zFD) - length(jFD);
    disp([' Removed ',num2str(rFD),' False Det']);
    zFD = jFD;
    if (isempty(zFD))
        zFD(1,1) = 1;  % keep from blowing up later
    end
    save(fn2,'zFD');
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find edges (start and end times) of bouts or sessions
dt = diff(ct)*24*60*60; % time between detections
%                   convert from days to seconds
I = [];
I = find(dt>gt);  % find start of gaps
sb = [ct(1);ct(I+1)];   % start time of bout
eb = [ct(I);ct(end)];   % end time of bout
dd = ct(end)-ct(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb);      % duration of bout in days
%find bouts > 10 sec long
% bd10 = find(bd > 1 / (60*60*24)); % for Kogia 10 sec
%disp(['Number Bouts : ',num2str(length(bd))])
% limit the length of a bout
blim = ltsamax/24;       % 6 hr bout length limit in days
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
bd = (eb - sb);   %duration bout in days
disp(['Number Bouts : ',num2str(nb)])
% LTSA session file
lspn = detpn;
lsfn = [sdn,'_',spe,'_LTSA',itnum,'.mat'];
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
% Test for False Detections
ixfd = (1: c4fd : length(ct));  % selected to test for False Det
A6 = exist(fn6,'file');
if (A6 ~= 2)
    zTD = -1.*ones(nb,4);
    save(fn6,'zTD');    % create new FD
    disp([' Make new TD2 file']);
else
    load(fn6)
    if (length(zTD(:,1)) ~= nb)
        disp([' Problem with TD file:',fn6]);
        return
    end
end
% Spectra Plot parameters
fimin = 0;
fimax = srate/2 ; % in kHz 100 or 160 kHz
% set ltsa step size
df = 1000*fimax/(length(pwr{1,1})-1);
% for LTSA PLOT
f = [1000*fimin:df:1000*fimax];
Ptf = interp1(tffreq,tfuppc,f,'linear','extrap'); % add to LTSA vector
if (specploton == 1) % do this only if spec plot
    % check length of MSP
    smsp = size(MSP);
    smsp2 = smsp(2); % 2nd element is num fft points
    for ift = 1 : 1: smsp2
        fmsp(ift) = ((srate/2)/(smsp2-1))*ift - (srate/2)/(smsp2-1);
    end
    fi = find(fmsp > fimin & fmsp <= fimax);
    fimint = fi(1); fimaxt = fi(end);
    fl = find(fmsp > flow);
    flowt = fl(1);
    ft = fmsp(fi);
    % for the PP vs RMS plot
    if (tfselect > 0)
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
cc = [' '];  % avoids crash when first bout too short
k = kstart;
% loop over the number of bouts (sessions)
while (k <= nb)
    disp([' BEGIN SESSION: ',num2str(k)]);
    % load in FD, MD and TD each session incase these have been modified
    load(fn2);  % brings in zFD
    %     load(fn3);  % brings in zMD
    load(fn6); % brings in zTD
    if (length(zTD(1,:)) == 2)
        zTD = [zTD,-1.*ones(length(zTD),2)];
        save(fn6,'zTD');
    end
    % Make PP versus RMS plot for all clicks
    if (specploton == 1)
        figure(51)
        hold off;
        xmpp = cl;
        i = find(xmpp > thres); % limit to high values
        Ndet = length(i);
        xmpp = xmpp(i);
        xmsp0 = csp(i,:) + ones(Ndet,1) *  Ptfpp;
        [xmsp,im] = max(xmsp0(:,flowt:fimaxt)'); % maximum between flow-100kHz
        for imax = 1 : 1 : length(im)
            Pmax = Ptfpp(im(imax) + flowt-1);
            xmpp(imax) = xmpp(imax) - tf + Pmax;
        end
        plot(xmsp,xmpp,'co')
        hold on;  % keep figure(51) plot with hold on
        xmpp = []; xmsp = [];
        title(['Based on ',num2str(length(i)),' clicks']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find detections and false detections within this bout (session)
    J = []; Jfandm =[]; Jtrue = []; XFD = [];
    J = find(ct >= sb(k) & ct <= eb(k));
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
    %
    if ~isempty(J)
        t = ct(J);        % detetion times in this session
        disp([' Detection times:',num2str(length(t))]);
        if (~isempty(XFD))
            xt = ct(ixfd(XFD));  %times to test for False Detection
            %xl = cl(XFD);   %amplitude for test False Detection
            disp([' Test False Detection times:',num2str(zTD(k,1))]),
        else
            xt = [];
        end
        y = cl(J);        % received levels in this session
        nd = length(J);     % number of detections in this session
        % get false detection times that intersect with detection times
        K2 = [];
        ff2 = [0];
        tfd = [];
        if (~isempty(zFD))
            [Cfd,K2,kfd] = intersect(ct(J),zFD(:,1));
        end
        if ~isempty(K2)
            ff2 = 1;
            tfd = ct(J(K2));
            yfd = cl(J(K2));
            ndS2 = length(K2);
            if (specploton == 1)
                if (ndS2 > 1)
                    wfd = sum(csn(J(K2),:))/ndS2;
                    sfd = sum(csp(J(K2),:))/ndS2;
                else   % must be one
                    wfd = csn(J(K2),:)/ndS2;
                    sfd = csp(J(K2),:)/ndS2;
                end
            end
            disp([' False detections:',num2str(ndS2)])
        else
            ff2 = 0;
            disp(' No False Det')
            ndS2 = 0;
        end
        Jfandm = [J(K2)];
        Jtrue = setxor(J,Jfandm);
        true = ct(Jtrue);
        ndS = length(true);
        if (specploton == 1)
            if (ndS > 1)
                wtrue = sum(csn(Jtrue,:))/ndS;
                strue = sum(csp(Jtrue,:))/ndS;
            else
                if (ndS == 1)
                    wtrue = csn(Jtrue,:);
                    strue = csp(Jtrue,:);
                end
            end
        end
        disp([' True Detections: ',num2str(ndS)])
    else
        disp('Error: no detections between bout start and end')
        return
    end
    dt = diff(t)*24*60*60; % inter-detection interval (IDI) and convert from days to seconds
    if ff2
        dtfd = dt(K2(1:end-1));
    end
    %     if ff3
    %         dtmd = dt(K3(1:end-1));
    %         %          disp([' Mis-ID times diff in this is session:',num2str(length(dtmd))])
    %     end
    disp(['END SESSION: ',num2str(k),' Start: ',datestr(sb(k)),...
        ' End:',datestr(eb(k))])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the number detection per bin
    dur = t(end) - t(1);    % session duration
    % sort detections into time bins, get max RL for time bin
    binDur = 5;     % bin duration [minutes] use 5 for density est
    nbin = ceil(dur*24*60/binDur);
    % interDetection Interval threshold
    dtTH = 1 ;  % [seconds] 1
    bin = 1;
    kb = 1;
    RL = zeros(nbin,1);
    C = zeros(nbin,1);  CX = zeros(nbin,1);
    T = zeros(nbin,1);
    Ndt = zeros(nbin,1);
    mdt = zeros(nbin,1);
    while kb <= nd
        %     while kb <= nbin
        tv = datevec(t(kb));     % put time(kb)in vector format
        tbin = floor(tv(5)/binDur);  % define time bin of time(kb)
        t0 = datenum([tv(1:4) tbin*binDur 00]); % lower bound
        t1 = datenum([tv(1:4) (tbin+1)*binDur 00]); % upper bound
        I = [];  IX = [];   % index for times in time bin
        I = find(t >= t0 & t < t1);  % find times in bin
        IX = find(xt >= t0 & xt < t1);  % find test times in bin
        if ~isempty(I)
            C(bin) = length(I); % the number of detections in time bin
            CX(bin) = length(IX); % the number of test detect in time bin
            RL(bin) = max(y(I));    % max RL in time bin
            T(bin) = t0;        % time for time bin
            if I > 1
                L = [];
                L = find(dt(I-1) < dtTH);
                if ~isempty(L)
                    Ndt(bin) = length(L);
                    % number of InterDetection Interval under threshold
                    mdt(bin) = mean(dt(I(L)-1));
                end
            end
            kb = I(end) + 1;      % set loop index to next time
            bin = bin + 1;  % increment bin number
        else
            if t(kb) == 0
                kb = kb + 1;
                bin = bin+1;
            else
                disp('this should not be possible')
                return
            end
        end
    end
    % filter empty and low number bins
    minNdet = 1;        % minimum number of detections per bin
    KB = [];  KBX = [];
    KB = find(C >= minNdet);
    KBX = find(CX >= minNdet);
    if isempty(KB)
        disp(['No bins with at least ',num2str(minNdet),' detections'])
        binT = [0];
        binRL = [0];
        binC = [0];
        k = k + 1;  % go to next
        continue
    else
        binT = T(KB) + datenum([0 0 0 0 binDur/2 0]);
        binRL = RL(KB);
        binC = C(KB);
        binCX = CX(KBX);
    end
    if (strcmp(cc,'w') && (zTD(k,2) == 0));
        disp(['Session: ',num2str(k),' # Test Detect Bins: ',...
            num2str(length(binCX)),' but NO False']);
        zTD(k,3) = length(binCX); zTD(k,4) = 0;
        save(fn6,'zTD');
        k = k + 1;
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number detection per Spectral bin in LTSA
    % make a spectra in figure 50
    PT = pt{1,k};   % LTSA session time vector
    pwr1 = pwr{1,k};  % LTSA power vector
    nbinS = length(PT);
    if (nbinS == 0)
        disp('No LTSA for this Session');
        PT(1) = sb(k) ; PT(2) = eb(k); % start end times for plots
        pwr1(1:length(f)) = ones; % make uniform LTSA
    else
        pwr1 = pwr1((1000*fimin/df)+1:(1000*fimax/df)+1,:); % limit plot range LTSA
    end
    durS = PT(end) - PT(1);
    minNdet = 1;        % minimum number of detections per bin
    %meanspec    %plot template mean spectra
    if (specploton == 1)
        if (ndS > 0)  % average true click spec
            mnSPEC = min(strue(fimint+5:flowt));
            SPEC = strue - mnSPEC;  % make low freq part = 0
            mxSPEC = max(SPEC(flowt:fimaxt));
            SPEC = SPEC ./ mxSPEC;  % make high freq part = 1
            figure(50);
            plot(ft,SPEC(fimint:fimaxt),'Linewidth',4)
            hold on;
            figure(52);
            plot(wtrue + 2*min(wtrue));
            hold on;
        else
            disp(['No true with at least ',num2str(minNdet),' detections'])
        end
        if (ndS2 > 0)  % average false click spec
            mnSPEC2 = min(sfd(fimint+5:flowt));
            SPEC2 = sfd - mnSPEC2;  % make low freq part = 0
            mxSPEC2 = max(SPEC2(flowt:fimaxt));
            SPEC2 = SPEC2 ./ mxSPEC2;  % make high freq part = 1
            figure(50);
            plot(ft,SPEC2(fimint:fimaxt),'r','Linewidth',4)
            figure(52);
            plot(wfd + min(wfd),'r');
        end
        figure(50)
        xlabel('Frequency (kHz)');
        ylim([0 1]);
        hold off;
        figure(52)
        xlabel(' 1ms total');
        hold off;
        figure(51)
        xlabel(' dB RMS '); ylabel(' dB PP ');
        % Add to PP versus RMS Plot for this session
        figure(51)
        Ndetpp = length(Jtrue);
        xmpp = cl(Jtrue);
        xmsp0 = csp(Jtrue,:) + ones(Ndetpp,1) *  Ptfpp;
        [xmsp,im] = max(xmsp0(:,flowt:fimaxt)');  % maximum value between 70 - 100 kHz
        for imax = 1 : 1 : length(im)
            Pmax = Ptfpp(im(imax) + flowt-1);
            xmpp(imax) = xmpp(imax) - tf + Pmax;
        end
        plot(xmsp,xmpp,'o')
        if (ndS2 > 0)  % add false click to figure (51) Jfandm
            Ndetppf = length(Jfandm);
            xmppf = cl(Jfandm);
            xmspf0 = csp(Jfandm,:) + ones(Ndetppf,1) *  Ptfpp;
            [xmspf,imf] = max(xmspf0(:,71:101)');  % maximum value between 70 - 100 kHz
            for imax = 1 : 1 : length(imf)
                Pmax = Ptfpp(imf(imax) + 70);
                xmppf(imax) = xmppf(imax) - tf + Pmax;
            end
            plot(xmspf,xmppf,'ro')
        end
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots stuff now in figure(201)
    warning('off')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(201)
    subplot(3,1,1)  % top panel RL vs Time
    plot(t,y,'.')
    hold on
    plot(t,y,'b.')
    if ff2 > 0
        plot(tfd,yfd,'r.')
        %         disp([' false det plotted:',num2str(length(tfd))])
    end
    hold off
    axis([PT(1) PT(end) rll rlh])
    datetick('x',15,'keeplimits')
    grid on
    tstr(1) = {fn};
    tstr(2) = {['Session: ',num2str(k),' Start Time ',...
        datestr(sb(k)),' Detect = ',num2str(nd)]};
    title(tstr);
    ylabel('RL [dB re 1\muPa]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,2)  % middle panel LTSA
    c = (contrast/100) .* pwr1 + bright;
    image(PT,f,c)
    axis xy
    v2 = axis;
    axis([PT(1) PT(end) v2(3) v2(4)])
    ylabel('Frequency Hz')
    datetick('keeplimits')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,1,3);  % bottom panel Inter-Detection Interval
    %dl = 1;  % dt limit seconds
    % make two copies of dt points for brush
    tdt2 = [];
    dt2 = [];
    if (length(dt) > 0)
        for idt = 1 : length(dt)
            tdt2 = [tdt2; t(idt) ; t(idt+1)];
            dt2 = [dt2 ; dt(idt) ; dt(idt) ];
        end
        [AX,H1,H2] = plotyy(tdt2,dt2,binT,binC,'plot','semilogy');
        set(H1,'Marker','.','MarkerFaceColor','b','LineStyle','none')
        set(H2,'Marker','o','MarkerFaceColor','c','LineStyle','none',...
            'Markersize',4.5)
        
        if ff2 > 0
            hold on
            plot(AX(1),tfd(2:end),dtfd,'.r') % no need to double FD
            %          since only the blue points are brush captured
            hold off
        end
        axis(AX(1),[PT(1) PT(end) 0 dl])
        axis(AX(2),[PT(1) PT(end) 1 100])
        datetick(AX(1),'x',15,'keeplimits')
        datetick(AX(2),'x',15,'keeplimits')
        Ytick = 0:dl/10:dl; % make 0.05 Kogia, 0.2 BW
        set(AX(1),'YTick',Ytick)
        
        Ytick2 = [.1 1 10 100 1000 10000];
        set(AX(2),'YTick',Ytick2)
        
        ylabel('Time between detections [s]')
        ylabel(AX(2),'Det/bin')
        xlabel('Time [GMT]')
        title('Inter-Detection Interval (IDI)')
        hold off
        grid(AX(1))
    else
        plot(0,0);
    end
    % end of plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % brush
    Bc = [];
    Bco = [];
    Bcs = [];
    pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get key stroke
    cc = get(201,'CurrentCharacter');
    if (strcmp(cc,'u') || ...
            strcmp(cc,'g') || strcmp(cc,'y') ||  strcmp(cc,'r') || ...
            strcmp(cc,'z') || strcmp(cc,'x') || strcmp(cc,'w'));
        disp(' Update Display') % Stay on same bout
    else
        if strcmp(cc,'s')
            dl = input(' Update IPI scale (sec):  '); % Set IPI scale
        else
            if strcmp(cc,'d')
                rll = input(' Update RL low (dB):  '); % Set RL low
                rlh = input(' Update RL high (dB):  '); % Set RL high
            else
                if strcmp(cc,'a')
                    contrast = input(' Update Contrast:  '); % Set RL high
                    bright = input(' Update Brightness:  '); % Set RL low
                else
                    if strcmp(cc,'b') % Move back one bout
                        if k ~= 1
                            k = k-1;
                        end
                    else
                        if strcmp(cc,'f') %assign as false
                            % False Detecctions
                            disp(['Number of False Detections = ',num2str(length(true))])
                            zFD = [zFD; true];   % cummulative False Detection matrix
                            save(fn2,'zFD')
                        else
                            if strcmp(cc,'t') %assign all as true
                                [C,iC] = setdiff(zFD(:,1),t);
                                disp(['Remaining False Detections = ',num2str(length(iC))])
                                zFD = zFD(iC,:);
                                save(fn2,'zFD')
                            else
                                if strcmp(cc,'j')% jump to seession
                                    prompt = 'Jump to Session: ';
                                    kjump = input(prompt);
                                    if (kjump > 0 && kjump < nb)
                                        k = kjump;
                                    end
                                else
                                    k = k+1;  % move forward one bout
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if (strcmp(cc,'x') || strcmp(cc,'z') ); % test click for random False Detect
        if (~isempty(XFD))
            zTD(k,2) = 0;
            for inxfd = 1 : zTD(k,1)
                figure(201)
                subplot(3,1,1)  % top panel RL vs Time
                hold on
                plot(ct(ixfd(XFD(inxfd))),cl(ixfd(XFD(inxfd))),...
                    'ro','MarkerSize',10);
                hold off;
                disp(['Showing #: ',num2str(inxfd),' click']);
                if (specploton == 1)
                    figure(50)  % add click to spec plot in BLACK
                    %meanspec;
                    plot(ft,SPEC(fimint:fimaxt),'Linewidth',4);
                    hold on;
                    tmnSPEC = min(csp(ixfd(XFD(inxfd)),fimint+5:flowt));
                    % make low freq part = 0
                    tSPEC = csp(ixfd(XFD(inxfd)),:) - tmnSPEC;
                    tmxSPEC = max(tSPEC(flowt:fimaxt));
                    tSPEC = tSPEC ./ tmxSPEC;  % make high freq part = 1
                    plot(ft,tSPEC(fimint:fimaxt),'k','Linewidth',4);
                    hold off;
                    %
                    figure(52) % add click to waveform plot in BLACK
                    plot(csn(ixfd(XFD(inxfd)),:),'k');
                    hold off;
                    pause
                    cc = get(52,'CurrentCharacter');
                    if (strcmp(cc,'z'))
                        zTD(k,2) = zTD(k,2) + 1;
                        zFD = [zFD; ct(ixfd(XFD(inxfd)))]; % add to FD
                    end
                end
            end
            disp([' Tested: ',num2str(zTD(k,1)),' False: ',...
                num2str(zTD(k,2))]);
            save(fn6,'zTD');
            save(fn2,'zFD')
        end
        k = k+1;
    end
    % Test 5 min window
    if (strcmp(cc,'w') && (zTD(k,2) > 0) ) ;  % test clicks in 5 min window
        % Plot all test clicks in session
        figure(201)
        subplot(3,1,1)  % top panel RL vs Time
        hold on
        for inxfd = 1 : zTD(k,1)
            plot(ct(ixfd(XFD(inxfd))),cl(ixfd(XFD(inxfd))),...
                'ro','MarkerSize',10);
        end
        hold off
        prompt = (['SESSION:',num2str(k),' #Test Detect Bins: ',...
            num2str(length(binCX)),' #False Bins: ']);
        pzTD = input(prompt);
        if (pzTD >= 0 && pzTD <= length(binCX) )
            zTD(k,4) = pzTD;
        else
            disp(' Number False Bins Out of Bounds')
            continue
        end
        zTD(k,3) = length(binCX);
        save(fn6,'zTD');
        k = k + 1;
        continue
    end
    % get brushed data from Figure(201) JAH 2-22-14
    hBrushLine = findall(gca,'tag','Brushing');
    brushData = get(hBrushLine, {'Xdata'});
    bsize = size(brushData);
    brushColor = get(hBrushLine, {'Color'});
    brushIdx = [];
    brushIdxo = [];
    brushIdxs = [];
    if ~isempty(brushData)
        brushIdx = ~isnan(brushData{1,1});  % get blue data
        %Check for Other Brush data
        if max(brushIdx) > 0
            % Put brush capture into Bc matrix
            Bc(:,1) = brushData{1,1}(brushIdx);
            if ((brushColor{1,1}(1,1) == 1.0 && brushColor{1,1}(1,2) == 0 ...
                    && brushColor{1,1}(1,3) == 0) || strcmp(cc,'r'));
                % Red paintbrush = False Detecctions
                disp(['Number of False Detections = ',num2str(length(Bc))])
                zFD = [zFD; Bc];   % cummulative False Detection matrix
                save(fn2,'zFD')
            else
                if ((brushColor{1,1}(1,1) == 1.0 &&brushColor{1,1}(1,2) == 1.0 ...
                        && brushColor{1,1}(1,3) == 0) || strcmp(cc,'y'));
                    % yellow paintbrush = give time of detection
                    disp(['       Year              Month          Day           Hour', ...
                        '          Min          Sec']);
                    disp(['Datevec ',num2str(datevec(Bc(1)))]);
                    yell = find(ct(J) == Bc(1));
                    if (specploton == 1)
                        figure(52) % add click to waveform plot in BLACK
                        hold off;
                        plot(wtrue + 2*min(wtrue),'c');
                        hold on;
                        plot(csn(J(yell),:),'k');
                        figure(50)  % add click to spec plot in BLACK
                        hold off
                        %meanspec;
                        plot(ft,SPEC(fimint:fimaxt),'Linewidth',4);
                        hold on
                        tmnSPEC = min(csp(J(yell),fimint+5:flowt));
                        % make low freq part = 0
                        tSPEC = csp(J(yell),:) - tmnSPEC;
                        tmxSPEC = max(tSPEC(flowt:fimaxt));
                        tSPEC = tSPEC ./ tmxSPEC;  % make high freq part = 1
                        plot(ft,tSPEC(fimint:fimaxt),'k','Linewidth',4);
                    end
                else
                    C =[];C2 = []; iC = []; iC2 = [];
                    disp(['Number of Detections Selected = ',num2str(length(Bc))])
                    if exist('zFD')
                        [C,iC] = setdiff(zFD(:,1),Bc(:,1));
                        disp(['Remaining Number of False Detections = ',num2str(length(iC))])
                        zFD = zFD(iC,:);
                        save(fn2,'zFD')
                    end
                end
            end
        end
    end
    Bc = [];
    if (specploton == 1)
        % get brush data for figure(52)
        figure(51)
        hBrush51 = findall(gca,'tag','Brushing');
        brushDataX = get(hBrush51, {'Xdata'});
        brushDataY = get(hBrush51, {'Ydata'});
        bsize51 = size(brushDataX);
        brushColor51 = get(hBrush51, {'Color'});
        brushId51x = [];  Cx = []; ixmpp = []; ibcy = [];
        brushId51y = [];  Bcx = [];  Bcy = [];
        if ~isempty(brushDataX)
            brushId51x = ~isnan(brushDataX{1,1});  % get dark blue data x
            brushId51y = ~isnan(brushDataY{1,1});  % get dark blue data y
            if max(brushId51y) > 0
                % Put brush capture into Bc matrix
                Bcx(:,1) = brushDataX{1,1}(brushId51x);
                Bcy(:,1) = brushDataY{1,1}(brushId51y);
            end
            [Cx, ixmpp, ibcy] = intersect(xmpp, Bcy);
            Bc = ct(Jtrue(ixmpp));  % assume false
            disp(['Number of False Detections = ',num2str(length(Bc))])
            zFD = [zFD; Bc];   % cummulative False Detection matrix
            save(fn2,'zFD')
        end
    end
    % don't end if you used paintbrush on last record
    if (k == nb && ~isempty(Bc))
        k = k-1;
        disp(' Last Record')
    end
    Bc = [];
    clear Bco;
    clear Bcs;
    clear hBrushLine; clear hBrush51;
    clear brushData; clear brushDataX; clear brushDataY;
    clear brushColor; clear brushColor51; clear bsize51;
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
