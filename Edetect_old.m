function Edetect_old
% Edetect.m
%
% function runs simple energy detector on xwav data.
% could be modified for typical wav files, but is xwav-centric
%
% user selects disk drives with xwav files/folders. The function
% reads all the xwav rawfile headers and finds rawfile times between user
% defined datenum times (ta and tb) from multiple events.
% Another function, Edetector (defined after this function) is called
% and runs a simple energy threshold detector from the time series composed
% of the multiple raw files between ta and tb.
% Based on energy detectors for clicks and sonar.
%
% User input: detector parameter xls pfile and event time period xls file
%
% Output parameters are saved to a binary matfile as
% detection times and peak amplitude levels
% 140919 adds addition outputs: snippet time series and spectra (from old Sdetect)
%150806 reads the transfer much and uses a dB threshold in the parameter
%file
%
%
% <warning: another suite of 'Edetect' code exist for LTSA sonar processing>
%
% initial build 130306 smw
%

clear all
clc

%
tic     % start execution clock
disp('Edetect_150806')
disp(' ')

% user interface get file via dialog box
% Parameter file
% [fname, pname] = uigetfile('*.txt','Load Parameter File');
disp('Load Parameter File');
[fname, pname] = uigetfile('*.xls;*.xlsx','Load Parameter File');
pfile = fullfile(pname, fname);
if strcmp(num2str(fname),'0')
    disp('Canceled Load Parameter File')
    return
else % give user some feedback
    disp(['Parameter File: ',pfile])
    disp(' ')
end
[num,txt] = xlsread(pfile);
fs = num(1); fa = num(2); fb = num(3); thDB = num(4); tffreq = num(5); ...
    lo = num(6); pl = num(7); un = num(8);
% fs = sample rate;
% fa fb = band pass filter start end frequency [Hz]
% thDB = 0-peak threshold [dB counts] for timeseries detections
% tffreq  = frequency from Transfer function to be applied to data
% lo = lock out time between consecutive detections [milliseconds]
% pl = pulse length to measure max ampliude and time [milliseconds]
% un = 1 for unfiltered data, 0 for no unfiltered

% user interface get file via dialog box
% Transfer Function file
% [fname, pname] = uigetfile('*.txt','Load Parameter File');
disp('Load Transfer Function File');
[tfFile, tfDir] = uigetfile('*.tf','Load Transfer Function File');
tfile = fullfile( tfDir, tfFile);
if strcmp(num2str(tfFile),'0')
    disp('Canceled Load Transfer Function File')
    return
else % give user some feedback
    disp(['Transfer Function: ',tfile])
    disp(' ')
end

%load transfer function
%N=512;
%fs = 200000;    %sample rate - should not be hardwired
% [tfFreq,tfPower] = textread(tfile,'%f %f'); %load transfer function
tfin = importdata(tfile);
tfFreq = tfin(:,1) ;  tfPower = tfin(:,2) ;
F=1:1:fs/2;
Ptf = interp1(tfFreq,tfPower,F,'linear','extrap');
% PtfN = downsample(Ptf,ceil(fs/N));
% set threshold in counts 0-peak
th = 10^((thDB - 6 - Ptf(tffreq))/20);

% user interface get file via dialog box
% Event Time file
% [fname, pname] = uigetfile([pname,'*.txt'],'Load Event Time File');
disp('Load event time file');
[fname, pname] = uigetfile([pname,'*.xls;*.xlsx'],'Load Event Time File');
tfile = fullfile(pname, fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save command line text to diary text file
[pathstr, name, ext] = fileparts(tfile);
dname = fullfile(pathstr,[name,'_Edetect_output.txt']);
diary(dname)
disp(dname)
disp(' ')
disp('Edetect_151202')
disp(' ')
disp(['Parameter File: ',pfile])
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(num2str(fname),'0')
    disp('Canceled Load Event Time File')
    return
else % give user some feedback
    disp(['Event Time File: ',tfile])
    disp(' ')
end
num = []; txt = [];
[num,txt] = xlsread(tfile);
% drl = txt(2:end,1);
% fd = txt(2:end,2);
toff = datenum([1900 0 -1 0 0 0]);   % xls is relative to 01Jan1900
dn1 = num(:,1) + toff; dn2 = num(:,2) + toff;
% dn1 = num(:,5) + toff; dn2 = num(:,6) + toff;   % for raw xls from Karli
% N = length(d1);
N = length(dn1);

disp(['Number of Events to run Edetector : ',num2str(N)])
disp(' ')

% [fd,d1,t1,d2,t2] = textread(tfile,'%s %s %s %s %s');
% tfile = path/filename for parameter file
% fd = path for xwav files to run on detector
% d1 t1 = date time start of bout 'mm/dd/yy' 'HH:MM:SS'
% d2 t2 = date time end of bout   'mm/dd/yy' 'HH:MM:SS'

% get drive letter of XWAV folders
drl = 1;
k = 1;
while ~strcmp(num2str(drl),'0')
    disp('Select Drive with Xwaves');
    drl = uigetdir('C:\','Select Drive with XWAV folders');
    dr{k} = drl;   % save disk names in cell array
    k = k + 1;
end
ndr = k - 2;    % number of disks

% load dnumstart and associated fullfilenames

% find folder names based on drive names provided
fk = 1; % folder counter
for k = 1: ndr          % loop over disks
    d = dir(dr{k});     % get directory list into structure d
    ndl = length(d);    % number of files and folders
    name = char(d.name);      % get name field of d structure
    for m = 1:ndl       % loop over file and folder names
        fldn = fullfile(dr{k},name(m,:));
        if isdir(fldn) && isempty(regexpi(fldn,'df'))
            fldrName{fk} = fldn;        % save folder names
            fk = fk + 1;
        end
    end
end
nf = fk - 1;
% char(fldrName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get dnumstart of all raw files in directories
NK = 0;  % counter for total number of raw files
fname = [];
for fk = 1 : nf      % loop over folders
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read xwav files names in directory
    d = [];
    d = dir(fullfile(deblank(char(fldrName{fk})),'*.x.wav'));    % xwav files
    
    fn = char(d.name);      % file names in directory
    fnsz = size(fn);        % number of data files in directory
    nfiles = fnsz(1);
    if nfiles < 1
        disp(['Error - No data files in this directory: ',fldrName{fk}])
        disp('Go to next directory')
        disp(' ')
        %     TT =0; PP = 0;
        continue
    else
        disp(' ')
        disp([num2str(nfiles),'  xwav data files in directory ',char(fldrName{fk})])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read times from each xwav header
    m = 0 ; % toal number of raw files of all xwavs in this directory
    for k = 1:nfiles
        fid = fopen(fullfile(deblank(char(fldrName{fk})),fn(k,:)),'r');
        fseek(fid,80,'bof');
        nrf = fread(fid,1,'uint16');         % Number of RawFiles in XWAV file (80 bytes from bof)
        fseek(fid,100,'bof');
        for r = 1:nrf                           % loop over the number of raw files in this xwav file
            m = m + 1;                                              % count total number of raw files
            %             rfileid(m) = r;                           % raw file id / number in this xwav file
            year(m) = fread(fid,1,'uchar');          % Year
            month(m) = fread(fid,1,'uchar');         % Month
            day(m) = fread(fid,1,'uchar');           % Day
            hour(m) = fread(fid,1,'uchar');          % Hour
            minute(m) = fread(fid,1,'uchar');        % Minute
            secs(m) = fread(fid,1,'uchar');          % Seconds
            ticks(m) = fread(fid,1,'uint16');        % Milliseconds
            %             byte_loc(m) = fread(fid,1,'uint32');     % Byte location in xwav file of RawFile start
            %             byte_length(m) = fread(fid,1,'uint32');    % Byte length of RawFile in xwav file
            %             write_length(m) = fread(fid,1,'uint32'); % # of blocks in RawFile length (default = 60000)
            %             sample_rate(m) = fread(fid,1,'uint32');  % sample rate of this RawFile
            %             gain(m) = fread(fid,1,'uint8');          % gain (1 = no change)
            %             padding = fread(fid,7,'uchar');    % Padding to make it 32 bytes...misc info can be added here
            %             fname(m,1:fnsz(2)) = fn(k,:);        % xwav file name for this raw file header
            
            NK = NK + 1;
            byte_loc(NK) = fread(fid,1,'uint32');     % Byte location in xwav file of RawFile start
            byte_length(NK) = fread(fid,1,'uint32');    % Byte length of RawFile in xwav file
            fseek(fid,16,'cof');
            fname{NK} = fullfile(deblank(char(fldrName{fk})),fn(k,:));        % xwav file name for this raw file header
            dnumStart(NK) = datenum([2000+year(m) month(m)...
                day(m) hour(m) minute(m) ...
                secs(m)+(ticks(m)/1000)]);
        end
        fclose(fid);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over events
%test N=2
% for k = 1:2
for k = 1:N
    %     ta = datenum([char(d1(k)),' ',char(t1(k))]);
    %     tb = datenum([char(d2(k)),' ',char(t2(k))]);
    ta = dn1(k);
    tb = dn2(k);
    %     dname = char(fd(k));
    %     dname = [char(drl(k)),char(fd(k))];
    disp(['Event ',num2str(k),'  ',datestr(ta,31),'  ',datestr(tb,31)])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find raw files during time period (clicking bout)
    I = [];
    I = find(ta <= dnumStart & tb > dnumStart);
    if ~isempty(I)
        if I(1) ~= 1
            I = [I(1)-1 I];     % need one before ta
            %         NRF = length(I);
            %         disp(['Number of rawfiles to evaluate ',num2str(NRF)])
        end
    else
        disp('Error: times not in directories or event shorter than 75s')
        disp(dr)
        TT{k} = ta; PP{k} = 0;
        SN{k} = ones(1,200+2);
        SP{k} = -10.*ones(1,101);
        USN{k} = ones(1,200+2);
        USP{k} = -10.*ones(1,101);
        NSP{k} = -10.*ones(1,101);
        continue
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [TT{k},PP{k},SN{k},USN{k},SP{k},USP{k},NSP{k}] = Edetector(I,fname,byte_loc,...
        byte_length,dnumStart,fa,fb,th,lo,pl,fs,tfFreq,tfPower);
    disp(['Number of Detections in Event = ',num2str(length(TT{k}))])
    disp(' ')
end

% move from cell to matrix
MTT = cell2mat(TT);
MPP = cell2mat(PP);
MSN = cell2mat(SN');
MSP = cell2mat(SP');
if (un == 1)
    MUSN = cell2mat(USN');
    MUSP = cell2mat(USP');
    MNSP = cell2mat(NSP');
end

% output file with Time and max RL
[pathstr, name, ext] = fileparts(tfile);
oname = fullfile(pathstr,[name,'_TPWS.mat']); % TimePeakWaveformSpectra
if (un == 1)
    save(oname,'MTT','MPP','MSN','MUSN','MSP','MUSP','MNSP')
else
    save(oname,'MTT','MPP','MSN','MSP')
end

disp(['Total Number of Detections for All Events = ',num2str(length(MTT))])
disp(' ')
t = toc;
disp(['Runtime : ',num2str(t),' seconds'])
diary off

function [TT,PP,SN,USN,SP,USP,NSP] = Edetector(I,fname,byte_loc,...
    byte_length,dnumStart,fa,fb,th,lo,pl,fs,tfFreq,tfPower)
%
% function runs simple energy detector on a group of individual raw files
% which define a single Event as per calling function Edetect
%
% Output:
% TT = pulse max amplitude detection time [datenum]
% PP = max amplitude peak-peak level[counts]
% 140919 add snippet time series and spectra for each detection
%
% Input:
% from calling function Edetect:
% I = indices for raw files beween Event start and end times
% fname = filenames 'vector' for rawfiles
% byte_loc = 'vector' of byte locations for rawfiles
% byte_length = 'vector' of number of bytes in rawfiles
% dnumStart = 'vecotr' of datenum start times for rawfiles
%
% from user input parameter file:
% fa = band pass filter low frequency corner [Hz]
% fb = band pass filter high frequency corner [Hz]
% th = threshold for detection [0-peak counts]
% lo = lock out time [milliseconds]
% pl = pulse length for max level measurement [milliseconds]
% fs = sample frequency
%
% stolen from detect1.m for sperm whale clicks - initially for beaked whale
% clicks on 4 chan tracking HARP, but also works for dolphin clicks.
% detector has a time lock out before next detection
%
% <warning: another suite of 'Edetect' code exist for LTSA sonar processing>
%
% initial build 130306 smw
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initialize
bps = 2;    % bytes per sample
nslo = fs * (lo/1000);  % number of samples to lock out before next detection
nspl = fs * (pl/1000);  % number of samples to use for max level measurement
[b,a] = ellip(4,0.1,40,[fa fb]*2/fs);   % elliptic filter

NRF = length(I);
disp(['Number of rawfiles to evaluate ',num2str(NRF)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over raw files within ta - tb event time period
%
TT = [];    % detection times
PP = [];    % pp [dB re counts] of detected pulse
SN = [];    % snippet timeseries of each detection
SP = [];    % spectra of each detection
USN = [];   % unfiltered snippet
USP = [];   % unfiltered spectra
NSP = [];   % noise spectra

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
durt = 0.001;  % seconds - snippet length
dur = durt * fs; % number samples for snippet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectra parameters
nfft = dur;
window = hanning(nfft);
overlap = 50;
noverlap = round((overlap/100)*nfft);


for k = 1:NRF
    R = I(k);   % raw file index within time period ta - tb
    %     disp(['Raw File ',num2str(R),'  ',fname(R,:),'   ',datestr(dnumStart(R),31)])
    disp(['Raw File ',num2str(R),'  ',fname{R}])
    %     fid = fopen(fullfile(dname,fname(R,:)),'r');
    fid = fopen(fname{R},'r');
    skip = byte_loc(R);
    fseek(fid,skip,-1);    % skip to correct rawfile
    spts = byte_length(R) / bps;
    ts = [];
    ts = fread(fid,spts,'int16');   % read one time slice ie one raw file
    fts = [];
    fts = filter(b,a,ts);
    %     fts = ts;
    
    % detect samples above threshold
    II = [];
    II = find(abs(fts) > th);   % sample number index of time series
    
    m = 1;      % counting index for II
    c = 1;      % counting index for detections
    J = [];     % sample index of detection of 1st sample > threshold
    T = [];     % detection time
    pp = [];    % max level within pulse length
    snip = [];  % matrix of timeseries snippets of each detection
    ufsnip = [];  % matrix of timeseries snippets of each detection
    noise = [];    % matrix of snip with click removed
    spec = [];  % matrix of spectras of each detection
    ufspec =[];
    nspec = [];
    
    % detect max amplitude of pulse and time
    while m <= length(II)
        K = [];     % I index after nslo samples from current (m) index
        K = find(II > II(m) + nslo,1,'first');  % find next trigger/detection after lock out
        if ~isempty(K)
            if II(m)+nspl <= spts
                [mx,M] = max(fts(II(m):II(m)+nspl));     % find max level within pl
                mn = min(fts(II(m):II(m)+nspl));               % find min level within pl
            else
                [mx,M] = max(fts(II(m):spts));     % find max level within pl
                mn = min(fts(II(m):spts));               % find min level within pl
            end
            pp(c) = mx-mn;     % peak to peak maximum
            % amplitude between the two triggered indices
            J(c) = II(m) + M - 1;  % index sample of max amplitude for this detection
            tds = J(c)/fs;
            T(c) = datenum([0 0 0 0 0 tds]) + dnumStart(R); % max amplitude detection datenum time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 140919 smw add timeseries snippet and spectra calc output
            
            if J(c) <= dur/2 +2  % put zeros on the front (ie beginning of file)
                dc = dur/2 +2 - J(c);
                if dc > 0
                    snip(c,:) = [zeros(dc,1);fts(1:J(c)+dur/2)]';
                    ufsnip(c,:) = [zeros(dc,1);ts(1:J(c)+dur/2)]';
                elseif dc == 0
                    snip(c,:) = fts(1:J(c)+dur/2);
                    ufsnip(c,:) = ts(1:J(c)+dur/2);
                end
            elseif J(c) + dur/2 >= spts % put zeros on the end (ie end of file)
                dc = spts - J(c) - dur/2;
                if dc == 0
                    snip(c,:) = [fts(J(c)-dur/2-1:spts)];
                    ufsnip(c,:) = [ts(J(c)-dur/2-1:spts)];
                else
                    snip(c,:) = [fts(J(c)-dur/2-1:spts);zeros(dc,1)]';
                    ufsnip(c,:) = [ts(J(c)-dur/2-1:spts);zeros(dc,1)]';
                end
            else
                snip(c,:) = fts(J(c)-dur/2-1:J(c)+dur/2);
                ufsnip(c,:) = ts(J(c)-dur/2-1:J(c)+dur/2);
            end
            [spec(c,:),f] = pwelch(snip(c,:),window,noverlap,nfft,fs);
            [ufspec(c,:),f] = pwelch(ufsnip(c,:),window,noverlap,nfft,fs);
            % make noise spectra
            noise = ufsnip;
            noise(c,150:200) = noise(c,1:51); % remove click 
            [nspec(c,:),f] = pwelch(noise(c,:),window,noverlap,nfft,fs);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            c = c + 1;
            m = K;  % move to next detections
        else
            m = m+1;
        end
    end
    Ndet = c - 1;
    disp(['Number of Detections in Raw file = ',num2str(Ndet)])
    %
    if (Ndet > 0)
        ppdb = 20.*log10(pp);     % logify
        %     snipdb = 20.*log10(snip);
        % apply transfer function
        specdb = 10.*log10(spec);
        ufspecdb = 10.*log10(ufspec);
        nspecdb = 10.*log10(nspec);
        Ptfx = interp1(tfFreq,tfPower,f,'linear','extrap');
        %         for i = 1:Ndet
        specdb = specdb + ones(Ndet,1) *  Ptfx';
        ufspecdb = ufspecdb + ones(Ndet,1) *  Ptfx';
        nspecdb = nspecdb + ones(Ndet,1) *  Ptfx';
        %         end
        % group for output
        TT = [TT T];
        PP = [PP ppdb];
        SN = [SN; snip];
        USN = [USN; ufsnip];
        SP = [SP; specdb];
        USP = [USP; ufspecdb];
        NSP = [NSP; nspecdb];
    end
    fclose(fid);    % close raw file
end
disp(' ')

