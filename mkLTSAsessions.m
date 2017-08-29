function mkLTSAsessions(varargin)
% mkLTSAsessions.m
% 2/21/15 version 1.1
% use individual click detections to define session/bout
% get and save ltsa pixel data for each session
% 140310 smw
% clear all
tic % start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set some parameters
gth = .5;    % gap time in hrs between sessions
gt = gth*60*60;    % gap time in sec

% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'filePrefix'
            filePrefix = varargin{n+1}; n=n+2;
        case 'detfn'
            detfn = varargin{n+1}; n=n+2;
        case 'sp'
            sp = varargin{n+1}; n=n+2;
        case 'lpn'
            lpn = varargin{n+1}; n=n+2;
        case 'sdir'
            sdir = varargin{n+1}; n=n+2;
        case 'srate'
            srate = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

fn = fullfile(sdir,detfn);
A1 = exist(fn,'file');
if A1 ~= 2
    disp(['Error: File Does Not Exist: ',fn])
    return
end

% ltsa file for GOM are named with GofMX and without underscore, and tf as
% well.
% E.g. GOM_DT_09 -> ltsa is GofMX_DT09
if findstr(detfn,'GOM')
    unscores = strfind(filePrefix,'_');
    filePrefix(unscores(2)) = '';
    filePrefix = regexprep(filePrefix,'GOM','GofMX');
end
%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:
p = sp_setting_defaults(sp,srate);

% user interface to get TF file
if (p.tfSelect > 0)
    disp('Load Transfer Function');
    if ~exist('tfName','var')% user interface to get TF file
        disp('Load Transfer Function File');
        [fname,pname] = uigetfile('I:\Harp_TF\*.tf','Load TF File');
        tffn = fullfile(pname,fname);
    else % or get it automatically from tf directory provided in settings
        stndeploy = strsplit(filePrefix,'_'); % get only station and deployment
        tffn = findTfFile(tfName,stndeploy); % get corresponding tf file
    end
    if strcmp(num2str(fname),'0')
        disp('Cancelled TF File');
        return
    else %give feedback
        disp(['TF File: ',tffn]);
    end
    fid = fopen(tffn);
    [A,~] = fscanf(fid,'%f %f',[2,inf]);
    tffreq = A(1,:);
    tfuppc = A(2,:);
    fclose(fid);
    
    tf = interp1(tffreq,tfuppc,p.tfSelect,'linear','extrap');
    disp(['TF @',num2str(p.tfSelect),' Hz =',num2str(tf)]);
else
    tf = 0;
    disp('No TF Applied ');
end

% LTSA session output file
lsfn = strrep(detfn,'TPWS','LTSA');
fn2 = fullfile(sdir,lsfn);

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
ib = find(cl >= p.threshRL);
disp([' Removed too low:',num2str(length(ia)-length(ib))]);
ct = ct(ib);
cl = cl(ib);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get ltsa file names for a specific site name and deployment number
d = dir(fullfile(lpn,[filePrefix,'*.ltsa']));
fnames = char(d.name);
nltsas = length(d);

% load up rawfile start times
doff = datenum([2000 0 0 0 0 0]);   % convert ltsa time to millenium time

% make the variables that hold the ltsa header info persistent, in case
% you are running itr_mkLTSA.m, this way you don't have to read the header
% info every time.
global sTime eTime rfTime

if isempty(sTime)
    if nltsas > 0
        sTime = zeros(nltsas,1); eTime = zeros(nltsas,1);
        disp('reading ltsa headers, please be patient ...')
        for k = 1:nltsas
            hdr = ioReadLTSAHeader(fullfile(lpn,fnames(k,:)));
            sTime(k) = hdr.ltsa.start.dnum + doff;  % start time of ltsa files
            eTime(k) = hdr.ltsa.end.dnum + doff;    % end time of ltsa files
            rfTime{k} = hdr.ltsa.dnumStart + doff; % all rawfiles times for all ltsas
        end
        disp('done reading ltsa headers')
    else
        disp(['No LTSAs found to match wildcard: ', [lpn,filePrefix,'*']])
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find edges (start and end times) of bouts or sessions
dt = diff(ct)*24*60*60; % time between detections
%                           convert from days to seconds
% I = [];
I = find(dt>gt);  % find start of gaps

if isempty(ct) % Catch in case there are no detections.
    return
end
sb = [ct(1);ct(I+1)];   % start time of bout
eb = [ct(I);ct(end)];   % end time of bout
% dd = ct(end)-ct(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb);      % duration of bout in days

% find bouts longer than the minimum
if ~isempty(p.minBout)
    bdI = find(bd > (p.minBout / (60*60*24)));
    bd = bd(bdI);
    sb = sb(bdI);
    eb = eb(bdI);
    nb = length(sb);        % number of bouts
end
% limit the length of a bout
blim = p.ltsaMax/24;       % 6 hr bout length limit
ib = 1;
while ib <= nb
    % disp([' ib = ',num2str(ib)]);
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

%Attempt to create a minimum bout duration
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1;
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
            L = find(rfTime{K} >= sb(k),1,'first');
        else
            L = find(rfTime{K} >= sb(k) & rfTime{K} <= eb(k));
        end
        if ~isempty(L)
            L = [L(1)-1,L]; % get rawfile from before sb(k)
            % grab the ltsa pwr matrix to plot
            hdr = ioReadLTSAHeader(fullfile(lpn,fnames(K,:))); % get some stuff we'll need
            nbin = length(L) * hdr.ltsa.nave(L(1));    % number of time bins to get
            fid = fopen(fullfile(hdr.ltsa.inpath,hdr.ltsa.infile),'r');
            % samples to skip over in ltsa file
            skip = hdr.ltsa.byteloc(L(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwr{k} = fread(fid,[hdr.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{K}(L(1));
            dt = datenum([0 0 0 0 0 5]);
            [ pt{k}, pwr{k} ] = padLTSAGaps(hdr, L,pwr{k});
            %   pt{k} = [t1:dt:t1 + (nbin-1)*dt]; % does not account for duty
            %   cycle
        else
            rfT = rfTime{K};
            disp('L is empty')
            pt{k} = [];
            pwr{k} = [];
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
        Ls = find(rfTime{Ks} >= sb(k));
        Le = [];
        Le = find(rfTime{Ke} <= eb(k));
        if ~isempty(Ls)
            Ls = [Ls(1)-1,Ls]; % get rawfile from before sb(k)
            % grab the ltsa pwr matrix to plot
            hdr = ioReadLTSAHeader(fullfile(lpn,fnames(Ks,:))); % get some stuff we'll need
            nbin = length(Ls) * hdr.ltsa.nave(Ls(1));    % number of time bins to get
            fid = fopen(fullfile(hdr.ltsa.inpath,hdr.ltsa.infile),'r');
            % samples to skip over in ltsa file
            skip = hdr.ltsa.byteloc(Ls(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwrLs = fread(fid,[hdr.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{Ks}(Ls(1));
            dt = datenum([0 0 0 0 0 5]);
            ptLs = [t1:dt:t1 + (nbin-1)*dt];
        end
        if ~isempty(Le)
            % grab the ltsa pwr matrix to plot
            hdr = ioReadLTSAHeader(fullfile(lpn,fnames(Ke,:))); % get some stuff we'll need
            nbin = length(Le) * hdr.ltsa.nave(Le(1));    % number of time bins to get
            fid = fopen(fullfile(hdr.ltsa.inpath,hdr.ltsa.infile),'r');
            % samples to skip over in ltsa file
            skip = hdr.ltsa.byteloc(Le(1));
            fseek(fid,skip,-1);    % skip over header + other data
            pwrLe = fread(fid,[hdr.ltsa.nf,nbin],'int8');   % read data
            fclose(fid);
            % make time vector
            t1 = rfTime{Ke}(Le(1));
            dt = datenum([0 0 0 0 0 5]);
            ptLe = t1:dt:t1 + (nbin-1)*dt;
        end
        
        if isempty(Ls) || isempty(Le)
            disp('Error: Ls or Le are empty ')
            k = k + 1;
            pwr{k} = [];
            pt{k} = [];
            continue
        else            % combine from end of ltsa with begin of next
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
% save(fn2,'pwr','pt')
save(fn2,'pwr','pt','-v7.3')
disp(['Done with file ',fn])
tc = toc;
disp(['Elasped Time : ',num2str(tc),' s'])

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                   Subroutines                           %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [ pt, pwr ] = padLTSAGaps(hdr,L,pwr)
%
% Takes LTSA raw file start times and returns a vector of start times for
% the individual spectral averages.  Accounts for duty cycle
%
% pt = vector containing start times of spectral averages
%
% hdr = Struct containing LTSA metadata
% L = vector of raw file indexes
%

pt = [];
mnum2secs = 24*60*60;
Y2K = datenum([ 2000 0 0 0 0 0 ]);
date_fmt = 'mm/dd/yy HH:MM:SS';

rfStarts = hdr.ltsa.dnumStart(L) + Y2K;
diffRF_s = round(diff(rfStarts)*mnum2secs);
rfDurs = unique(diffRF_s);
nbin = length(L) * hdr.ltsa.nave(L(1));
t1 = rfStarts(1);
dt_dnum = datenum([ 0 0 0 0 0 hdr.ltsa.tave ]);
fs = hdr.ltsa.fs;

rfdt = (hdr.ltsa.dnumEnd(1) - hdr.ltsa.dnumStart(1))*(24*60*60);
% if fs == 200e3
%     rfdt = 75; % 75 second raw file duration
% else
%     fprintf('Unsupported sample rate of %d\n', fs);
%     mfn = mfilename('fullpath');
%     fprintf('Add raw file duration [seconds] for sample rate to %s\n', mfn);
%     fprintf('Or ask JAH...definitely ask JAH\n');
% end


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
            numaves2pad = ceil(tdt/hdr.ltsa.tave);
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