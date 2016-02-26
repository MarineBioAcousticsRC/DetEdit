function mkLTSAsessions(~)
% mkLTSAsessions.m
% 2/21/15 version 1.1
% use individual click detections to define session/bout
% get and save ltsa pixel data for each session
% 140310 smw
clear all
tic % start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some parameters
gth = .5;    % gap time in hrs between sessions
gt = gth*60*60;    % gap time in sec
ltsamax = 6; % length of ltsa window
thres = 80;
% get user input and set up file names
stn = input('Enter Project Name (MC GC DT HH JAX SOCAL): ','s'); % site name
dpn = input('Enter Deployment number (01 02 11D ...): ','s'); % deployment number
itnum = input('Enter Iteration number (1 2 ...): ','s');
sdn = [stn,dpn];    % site name and deployment number
sp = input('Enter Species: Zc Me BWG Md Ko De ','s');
if (strcmp(sp,'Ko') || strcmp(sp,'k'))
    specchar = 'K'; %Simone abbreviation for species
    spe = 'Kogia';  tfselect = 80000; % freq used for transfer function
elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
    specchar = 'Z'; %Simone abbreviations for BW species
    spe = 'Cuviers';  tfselect = 40200; % freq used for transfer function
    thres = 121;
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    specchar = 'M'; %Simone abbreviations for BW species
    spe = 'Gervais'; tfselect = 40200; % freq used for transfer function
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    specchar = 'G'; %Simone abbreviations for BW species
    spe = 'BWG';  tfselect = 40200; % freq used for transfer function
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    specchar = 'D'; %Simone abbreviations for BW species
    spe = 'BW31';  tfselect = 40200; % freq used for transfer function
elseif (strcmp(sp,'De') || strcmp(sp,'de'))
    %specchar = 'D'; %Simone abbreviations for BW species
    spe = 'Delphin';  tfselect = 0; % freq used for transfer function
elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
    %specchar = 'D'; %Simone abbreviations for BW species
    spe = 'Porpoise';  tfselect = 0; % freq used for transfer function
    thres = 100;
elseif (strcmp(sp,'MFA') || strcmp(sp,'mfa'))
    spe = 'MFA';  tfselect = 4000; % freq used for transfer function
    ltsamax = .5; % length of ltsa window
    thres = 80;
elseif (strcmp(sp,'Dl') || strcmp(sp,'dl'))
    spe = 'Beluga'; tfselect = 45000; %freq used for tf
else
    disp(' Bad Species type')
    return
end
disp('Select Directory with LTSA');
lpn = uigetdir('I:\','Select Directory with LTSAs');
lpn = [lpn,'\'];
% lpn = ['I:\JAH\ltsa\',stn,'\'];  % ltsa pathname
% detection file
disp('Select Directory with Detections');
sdir = uigetdir('I:\','Select Directory with Detections');
%
% detpn = [sdir,['\',stn,'_',spe],'\'];
detpn = [sdir,'\'];
% detfn = [stn,dpn,'_',spe,'_TPWS2.mat'];
detfn = [stn,dpn,'_',spe,'_TPWS',itnum,'.mat'];
fn = fullfile(detpn,detfn);
A1 = exist(fn,'file');
if A1 ~= 2
    disp(['Error: File Does Not Exist: ',fn])
    return
end
% user interface to get TF file
if (tfselect > 0)
    disp('Load Transfer Function');
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
    
    tf = interp1(tffreq,tfuppc,tfselect,'linear','extrap');
    disp(['TF @',num2str(tfselect),' Hz =',num2str(tf)]);
else
    tf = 0;
    disp('No TF Applied ');
end
% LTSA session output file
lspn = detpn;
lsfn = [sdn,'_',spe,'_LTSA',itnum,'.mat'];
%lsfn = [sdn,'_',spe,'_LTSA.mat'];
fn2 = fullfile(lspn,lsfn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load detections
% MTT = click detection time, MPP = click RL level (dBpp counts)
% MSN = Nclicks x length snip, MSP = Nclicks x spectra
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
% apply tf and remove low amplitude
cl = cl + tf;
ib = find(cl >= thres);
disp([' Removed too low:',num2str(length(ia)-length(ib))]);
ct = ct(ib);
cl = cl(ib);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get ltsa file names for a specific site name and deployment number
%d = dir([lpn,'GofMX_',sdn,'*']);
if (strcmp(stn,'MC')|| strcmp(stn,'GC') || strcmp(stn,'DC') || ...
        strcmp(stn,'DT') || strcmp(stn,'HH'))
    d = dir([lpn,'GofMX_',sdn,'*']);
else
    d = dir([lpn,sdn,'*']);
end
fnames = char(d.name);
nltsas = length(d);
% load up rawfile start times
disp('reading ltsa headers, please be patient ...')
doff = datenum([2000 0 0 0 0 0]);   % convert ltsa time to millenium time
global PARAMS
sTime = zeros(nltsas,1); eTime = zeros(nltsas,1);
for k = 1:nltsas
    % for k = 1:2
    PARAMS.ltsa.inpath = lpn;
    PARAMS.ltsa.infile = fnames(k,:);
    read_ltsahead_GoM
    sTime(k) = PARAMS.ltsa.start.dnum + doff;  % start time of ltsa files
    eTime(k) = PARAMS.ltsa.end.dnum + doff;    % end time of ltsa files
    
    rfTime{k,:} = PARAMS.ltsa.dnumStart + doff; % all rawfiles times for all ltsas
end
% hack fix
if strcmp(sdn,'GC02')
    eTime(2) = sTime(3);
end
disp('done reading ltsa headers')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find edges (start and end times) of bouts or sessions
dt = diff(ct)*24*60*60; % time between detections
%                           convert from days to seconds
I = [];
I = find(dt>gt);  % find start of gaps
sb = [ct(1);ct(I+1)];   % start time of bout
eb = [ct(I);ct(end)];   % end time of bout
dd = ct(end)-ct(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb);      % duration of bout in days
%find bouts > 10 sec long
% bd10 = find(bd > 1 / (60*60*24)); % for Kogia 10 sec
% disp(['Number Bouts: ',num2str(length(bd))])
% limit the length of a bout
blim = ltsamax/24;       % 6 hr bout length limit
ib = 1;
while ib <= nb
    %disp([' ib = ',num2str(ib)]);
    bd = (eb - sb);   %duration bout in sec
    if (bd(ib) > blim)      % find long bouts
        nadd = ceil(bd(ib)/blim) - 1; % number of bouts to add
        for imove = nb : -1: (ib +1)
            sb(imove+nadd)= sb(imove);
            %disp(['imove-sb  ',num2str(imove)])
        end
        for iadd = 1 : 1: nadd
            sb(ib+iadd) = sb(ib) + blim*iadd;
            %disp(['iadd-sb',num2str(iadd)])
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
bd = (eb - sb);   %duration bout in sec
disp(['Number Bouts : ',num2str(nb)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
% for debugging only
% k = 200;fprintf('DEBUGGING VARIABLE IN PLACE, STARTING AT SESSIONS %d\n',k);
% loop over the number of bouts (sessions)
while (k <= nb)
    %      if eb(k) - sb(k) < 1 / (60*60*24)
    %         disp('Session less than 1s, so skip it')
    %         k = k + 1;
    %         continue
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find which ltsa to use and get pwr and pt vectors
    K = [];
    K = find(sTime <= sb(k) & eTime >= eb(k));
    % find which rawfiles to plot ltsa
    if ~isempty(K) && length(K) == 1
        L = [];
        if eb(k) -sb(k) < 75/ (60*60*24)
            L = find(rfTime{K,:} >= sb(k),1,'first');
        else
            L = find(rfTime{K,:} >= sb(k) & rfTime{K,:} <= eb(k));
        end
        if ~isempty(L)
            L = [L(1)-1,L]; % get rawfile from before sb(k)
            % grab the ltsa pwr matrix to plot
            PARAMS.ltsa.inpath = lpn;
            PARAMS.ltsa.infile = fnames(K,:);
            read_ltsahead_GoM % get some stuff we'll need
            nbin = length(L) * PARAMS.ltsa.nave(L(1));    % number of time bins to get
            fid = fopen([PARAMS.ltsa.inpath,PARAMS.ltsa.infile],'r');
            % samples to skip over in ltsa file
            skip = PARAMS.ltsa.byteloc(L(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwr{k} = fread(fid,[PARAMS.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{K}(L(1));
            dt = datenum([0 0 0 0 0 5]);
            [ pt{k}, pwr{k} ] = padLTSAGaps(PARAMS, L,pwr{k});
            %   pt{k} = [t1:dt:t1 + (nbin-1)*dt]; % does not account for duty
            %   cycle
        else
            rfT = rfTime{K,:};
            disp('L is empty')
            %             disp(datestr(rfT))
            disp(['bout start time is ',datestr(sb(k))])
            disp(['bout end time is ',datestr(eb(k))])
        end
    elseif isempty(K)   % use the end of one and the beginning of next ltsa
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
        Ls = find(rfTime{Ks,:} >= sb(k));
        Le = [];
        Le = find(rfTime{Ke,:} <= eb(k));
        if ~isempty(Ls)
            Ls = [Ls(1)-1,Ls]; % get rawfile from before sb(k)
            % grab the ltsa pwr matrix to plot
            PARAMS.ltsa.inpath = lpn;
            PARAMS.ltsa.infile = fnames(Ks,:);
            read_ltsahead_GoM % get some stuff we'll need
            nbin = length(Ls) * PARAMS.ltsa.nave(Ls(1));    % number of time bins to get
            fid = fopen([PARAMS.ltsa.inpath,PARAMS.ltsa.infile],'r');
            % samples to skip over in ltsa file
            skip = PARAMS.ltsa.byteloc(Ls(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwrLs = fread(fid,[PARAMS.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{Ks}(Ls(1));
            dt = datenum([0 0 0 0 0 5]);
            ptLs = [t1:dt:t1 + (nbin-1)*dt];
        end
        if ~isempty(Le)
            % grab the ltsa pwr matrix to plot
            PARAMS.ltsa.inpath = lpn;
            PARAMS.ltsa.infile = fnames(Ke,:);
            read_ltsahead_GoM % get some stuff we'll need
            nbin = length(Le) * PARAMS.ltsa.nave(Le(1));    % number of time bins to get
            fid = fopen([PARAMS.ltsa.inpath,PARAMS.ltsa.infile],'r');
            % samples to skip over in ltsa file
            skip = PARAMS.ltsa.byteloc(Le(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwrLe = fread(fid,[PARAMS.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{Ke}(Le(1));
            dt = datenum([0 0 0 0 0 5]);
            ptLe = [t1:dt:t1 + (nbin-1)*dt];
        end
        
        if isempty(Ls) || isempty(Le)
            disp('Error: Ls or Le are empty ')
            k = k + 1;
            continue
        else            % combine from end of ltsa with begin of next
            pwr{k} = [pwrLs pwrLe];
            pt{k} = [ptLs ptLe];
        end
    else
        disp(['K = ',num2str(K)])
        disp(['bout start time is ',datestr(sb(k))])
        disp(['bout end time is ',datestr(eb(k))])
    end
    
    disp(['Session: ',num2str(k),'  Start: ',datestr(sb(k)),'  End:',datestr(eb(k)),...
        '   Duration: ',num2str(24*60*60*bd(k)),' sec'])
    
    k = k+1;
end
% save(fn2,'pwr','pt')
save(fn2,'pwr','pt','-v7.3')
disp(['Done with file ',fn])
tc = toc;
disp(['Elasped Time : ',num2str(tc),' s'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                   Subroutines                           %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ pt, pwr ] = padLTSAGaps(PARAMS,L,pwr)
%
% Takes LTSA raw file start times and returns a vector of start times for
% the individual spectral averages.  Accounts for duty cycle
%
% pt = vector containing start times of spectral averages
%
% PARAMS = Struct containing LTSA metadata
% L = vector of raw file indexes
%

pt = [];
mnum2secs = 24*60*60;
Y2K = datenum([ 2000 0 0 0 0 0 ]);
date_fmt = 'mm/dd/yy HH:MM:SS';

rfStarts = PARAMS.ltsa.dnumStart(L) + Y2K;
diffRF_s = round(diff(rfStarts)*mnum2secs);
rfDurs = unique(diffRF_s);
nbin = length(L) * PARAMS.ltsa.nave(L(1));
t1 = rfStarts(1);
dt_dnum = datenum([ 0 0 0 0 0 PARAMS.ltsa.tave ]);
fs = PARAMS.ltsa.fs;

if fs == 200e3
    rfdt = 75; % 75 second raw file duration
else
    fprintf('Unsupported sample rate of %d\n', fs);
    mfn = mfilename('fullpath');
    fprintf('Add raw file duration [seconds] for sample rate to %s\n', mfn);
    fprintf('Or ask JAH...definitely ask JAH\n');
end


nfreq = size(pwr,1);
totalaves = size(pwr,2);
gapsidx = find(diffRF_s~=rfdt);
rfdt_dnum = datenum([ 0 0 0 0 0 rfdt-5 ]);
if size(rfDurs, 2) > 1
    fprintf('Duty cycle or gap encountered!\n')
    fprintf('\tBout starting at %s\n', datestr(rfStarts(1),date_fmt));
    pt = t1:dt_dnum:t1+rfdt_dnum; % make first raw file's part of pt;
    rf = 2;
    while rf <= length(L)
        if ismember(rf,gapsidx+1) % if there's a duty cycle or gap, pad pwr for later
            tdt = round((rfStarts(rf) - rfStarts(rf-1))*mnum2secs-rfdt);
            numaves2pad = ceil(tdt/PARAMS.ltsa.tave);
            pwrpad = ones(nfreq,numaves2pad).*-128;
            gaveidx = length(pt); % last average before gap
            pwr = [ pwr(:,1:gaveidx), pwrpad, pwr(:,gaveidx+1:end)];
            padavetimes = pt(end)+dt_dnum:dt_dnum:rfStarts(rf)-dt_dnum;
            pt = [ pt, padavetimes ];
        end
        
        tpt = rfStarts(rf):dt_dnum:rfStarts(rf)+rfdt_dnum;
        pt = [ pt, tpt ];
        rf = rf+1;
    end
else
    pt = [t1:dt_dnum:t1 + (nbin-1)*dt_dnum];
end

end