function Edetect(varargin)
% Edetect.m
%
% function runs a simple energy detector on audio files (.wav or .x.wav).
%
% user selects disk drives with audio files/folders. The function
% reads all the audio file headers and finds files times between user
% defined datenum times (ta and tb) from multiple events.
% Another function, Edetector (defined after this function) is called
% and runs a simple energy threshold detector from the time series composed
% of the multiple audio files between ta and tb.
% Based on energy detectors for clicks and sonar.
%
% User input: detector parameter xls pfile and event time period xls file
% Input can be done through GUI tools, if not provided in function call.
% All arguments are optional.
% Example function call:
%   Edetect_wav('paramFile','G:\EndToEnd_Example_Test_Set\edetect_test_params.xlsx',...
%         'tfFile','E:\Code\TF_files\HF631_140122\HF631_140122_sig1_invSensit.tf',...
%         'timeFile','G:\GofMXArraySpRecs\Pc\MOTU\edetect_test.xlsx',...
%         'audioDir',{'G:\GofMXArraySpRecs\Pc\MOTU';'G:\GofMXArraySpRecs\Gmsp\MOTU'},...
%         'fileType',1,...
%         'channel',1,...
%         'outputFileName','siteA_2012')
%   paramFile: Spreadsheet containing detector parameters
%   tfFile: Transfer function
%   timeFile: Spreadsheet containing bout times
%   audioDir: Folder name, or cell array of multiple folder names to search
%       for wav files. NOTE: Search includes this folder specified and any
%       subfolders (eg: 'E:\myFolder' and 'E:\myFolder\mySubFolder').
%   fileType: Audio file type: 1 - wav, 2 - pseudo-wav (x.wav). Defaults to 2.
%   channel: Which wav file channel to detect on. Defaults to 1.
%   outputFileName: Output File Name


% Output parameters are saved to a binary matfile as
% detection times and peak amplitude levels

% 140919 adds addition outputs: snippet time series and spectra (from old Sdetect)
% 150806 reads the transfer much and uses a dB threshold in the parameter
% file
%
%
% <warning: another suite of 'Edetect' code exist for LTSA sonar processing>
%
% initial build 130306 smw
% wav adaptation 17-020 kef

clc

tic     % start execution clock
fprintf('Edetect\n\n')


% check for input arguments
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'paramFile'
            paramFile = varargin{n+1}; n=n+2;
        case 'tfFile'
            tfFile = varargin{n+1}; n=n+2;
        case 'timeFile'
            timeFile = varargin{n+1}; n=n+2;
        case 'audioDir'
            audioDir = varargin{n+1}; n=n+2;
            if ~iscell(audioDir)
                % expects a cell array. If no cell array, put the string
                % into a single cell array.
                audioDirVec = audioDir;
                audioDir = cell(1);
                audioDir{1} = audioDirVec;
            end
        case 'fileType'
            ftype = varargin{n+1}; n=n+2;
        case 'channel'
            channel = varargin{n+1}; n=n+2;
        case 'outputFileName'
            outputFileName = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

% Import detection parameters from spreadsheet
if ~exist('paramFile','var')
    % create user interface get file via dialog box
    [paramFileName, paramFilePath] = uigetfile('*.xls;*.xlsx','Load Parameter File');
    paramFile = fullfile(paramFilePath, paramFileName);
    if strcmp(num2str(paramFilePath),'0')
        disp('Canceled Load Parameter File')
        return
    end
end

fprintf('Parameter File: %s\n',paramFile)
disp('Loading Parameter File');
[num,~] = xlsread(paramFile);
fs = num(1); fa = num(2); fb = num(3); thDB = num(4); tfAdjustFreq = num(5); ...
    lo = num(6); pl = num(7); un = num(8);
% Definitions:
% fs = sample rate;
% fa fb = band pass filter start end frequency [Hz]
% thDB = 0-peak threshold [dB counts] for timeseries detections
% tffreq  = frequency from Transfer function to be applied to data
% lo = lock out time between consecutive detections [milliseconds]
% pl = pulse length to measure max ampliude and time [milliseconds]
% un = 1 for unfiltered data, 0 for no unfiltered


% Import transfer function
if ~exist('tfFile','var')
    % create user interface get transfer function file via dialog box
    [tfFileName, tfDir] = uigetfile('*.tf','Load Transfer Function File');
    tfFile = fullfile(tfDir, tfFileName);
    if strcmp(num2str(tfFileName),'0')
        disp('Canceled Load Transfer Function File')
        return
    end
end

fprintf('Transfer Function: %s\n',tfFile)
fprintf('Loading Transfer Function File\n\n')
tfInput = importdata(tfFile);
tfFreq = tfInput(:,1) ;
tfPower = tfInput(:,2) ;
F = 1:1:fs/2;
Ptf = interp1(tfFreq,tfPower,F,'linear','extrap');
% PtfN = downsample(Ptf,ceil(fs/N));
% set threshold in counts 0-peak
countThresh = 10^((thDB - 6 - Ptf(tfAdjustFreq))/20);


% Get event times from spreadsheet
if ~exist('timeFile','var')
    % create user interface get file via dialog box
    [timeFileName, timeFileDir] = uigetfile('*.xls;*.xlsx','Load Event Time File');
    timeFile = fullfile(timeFileDir, timeFileName);
    if strcmp(num2str(timeFile),'0')
        disp('Canceled Load Time Event File')
        return
    end
end
fprintf('Event Time File: %s\n',timeFile)
fprintf('Load event time file\n\n');
[dates2Detect,~] = xlsread(timeFile);

if ~exist('fileType','var')
    fprintf('No file type specified. Reading pseudo-wav (x.wav) files.\n\n')
    ftype = 2;
end

if ~exist('channel','var')
    fprintf('No channel specified. Detecting on channel 1.\n\n')
    channel = 1;
end

if ~exist('outputFileName','var')
   [~, outputFileName, ~] = fileparts(timeFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save command line text to diary text file
[pathstr, ~, ~] = fileparts(timeFile);
dname = fullfile(pathstr,[outputFileName,'_Edetect_output.txt']);
diary(dname)
fprintf('%s\n\n',dname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% drl = txt(2:end,1);
% fd = txt(2:end,2);
timeOffset = datenum([1900 0 -1 0 0 0]);   % xls is relative to 01Jan1900
boutStartDNum = dates2Detect(:,1) + timeOffset;
boutEndDNum = dates2Detect(:,2) + timeOffset;
N = length(boutStartDNum);

fprintf('Number of Events to run Edetector : %.0f\n\n', N)

% [fd,d1,t1,d2,t2] = textread(tfile,'%s %s %s %s %s');
% tfile = path/filename for parameter file
% fd = path for xwav files to run on detector
% d1 t1 = date time start of bout 'mm/dd/yy' 'HH:MM:SS'
% d2 t2 = date time end of bout   'mm/dd/yy' 'HH:MM:SS'

% get drive letter of XWAV folders
if ~exist('audioDir','var')
    audioDir = 1;
    k = 1;
    while ~strcmp(num2str(audioDir),'0')
        disp('Select folder with audio files');
        audioDir = uigetdir('C:\','Select folder with wav folders');
        audioDirList{k} = audioDir;   % save folder names in cell array
        k = k + 1;
    end
    ndr = k- 2;    % number of folders
    
else
    k = length(audioDir);
    audioDirList = audioDir;
    ndr = k;
end

% load dnumstart and associated fullfilenames

% find folder names based on drive names provided
fk = 1; % folder counter
for k = 1: ndr          % loop over disks
    fldrName{fk} = audioDirList{k}; % add base directory to folder list
    fk = fk + 1;
    
    % Now check for subfolders that might contain .wav files
    subDirList = dir([audioDirList{k},'/']);     % get subdirectory list into structure d
    subDirList = subDirList(3:end); % get rid of the . and .. folders
    isDir = [subDirList.isdir];
    subDirListPrune = subDirList(isDir); % remove anythign that isn't a directory
    ndl = length(subDirListPrune);    % number of files and folders
    name = char(subDirListPrune.name);      % get name field of d structure
    for m = 1:ndl       % loop over file and folder names
        fldn = fullfile(audioDirList{k}, subDirListPrune(m).name);
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
DateRE = '\d+_\d+';
fullFileName = [];
for fk = 1 : nf      % loop over folders
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find audio file names in directory
    d = [];
    if ftype == 1
        d = dir(fullfile(deblank(char(fldrName{fk})),'*.wav'));    % wav files
    elseif ftype == 2
        d = dir(fullfile(deblank(char(fldrName{fk})),'*x.wav'));    % xwav files
    end
    
    fn = char(d.name);      % file names in directory
    fnsz = size(fn);        % number of data files in directory
    nfldrFiles = fnsz(1);
    if nfldrFiles < 1
        disp(['No data files in this directory: ',fldrName{fk}])
        disp('Go to next directory')
        disp(' ')
        %     TT =0; PP = 0;
        continue
    else
        disp(' ')
        disp([num2str(nfldrFiles),'  audio data files in directory ',char(fldrName{fk})])
        for w = 1:nfldrFiles % build full path for each file
            fullFileName = [fullFileName;fullfile(fldrName{fk},fn(w,:))];
        end
    end
end
nfiles = size(fullFileName,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over bouts from event time files
for b = 1:N
    % read times from each audio file name
    % loop over audio files and run detector
    for k = 1:nfiles
        inFileName = fullfile(fullFileName(k,:));
        hdrs(k) = ioReadXwavHeader(inFileName,DateRE,ftype);
        % check if this file overlaps with a time we are interested in
        thisFileStart = hdrs(k).start.dnum;
        thisFileEnd = hdrs(k).end.dnum;
        afterEndList = thisFileStart>=boutEndDNum(b);
        beforeStartList = thisFileEnd<=boutStartDNum(b);
        % if sum of these two vectors = 0 anywhere, it means the file
        % overlaps with a time we're interested in.
        overlapBouts = find(afterEndList+beforeStartList==0);
        
        if ~isempty(overlapBouts)
            [TT{k},PP{k},SN{k},USN{k},SP{k},USP{k},NSP{k}] = Edetector(hdrs(k),...
                inFileName,fa,fb,countThresh,lo,pl,tfFreq,tfPower,channel,ftype,...
                boutStartDNum(b),boutEndDNum(b));
            fprintf('Number of Detections in Event = %.0f\n\n',length(TT{k}))
            
        else
            fprintf('Skipping file %s\n',inFileName)
            fprintf('Outside bout ranges.\n\n')
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% move from cell to matrix
MTT = cell2mat(TT)';
MPP = cell2mat(PP)';
MSN = cell2mat(SN);
MSP = cell2mat(SP);
if (un == 1)
    MUSN = cell2mat(USN');
    MUSP = cell2mat(USP');
    MNSP = cell2mat(NSP');
end

% output file with Time and max RL
outputPath = fullfile(pathstr,[outputFileName,'_TPWS1.mat']); % TimePeakWaveformSpectra
if (un == 1)
    save(outputPath,'MTT','MPP','MSN','MUSN','MSP','MUSP','MNSP')
    disp(['TPWS file saved in ',outputPath])
else
    save(outputPath,'MTT','MPP','MSN','MSP')
    disp(['TPWS file saved in ',outputPath])
end

disp(['Total Number of Detections for All Events = ',num2str(length(MTT))])
disp(' ')
t = toc;
disp(['Runtime : ',num2str(t),' seconds'])
diary off

function [TT,PP,SN,USN,SP,USP,NSP] = Edetector(hdr,fname,...
    fa,fb,th,lo,pl,tfFreq,tfPower,channel,ftype,startBout,endBout)
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
fs = hdr.fs;
bps = 2;    % bytes per sample
nslo = hdr.fs * (lo/1000);  % number of samples to lock out before next detection
nspl = fs * (pl/1000);  % number of samples to use for max level measurement
[b,a] = ellip(4,0.1,40,[fa fb]*2/fs);   % elliptic filter

[startsSec,stopsSec] = dST_choose_segments(hdr,ftype);

% find chunks of data during bout
% convert bout times in seconds from start of audio file
startBoutSec = (startBout - hdr.start.dnum)*24*60*60;
endBoutSec = (endBout - hdr.start.dnum)*24*60*60;
ichunk = find(startBoutSec <= startsSec & endBoutSec > startsSec);

startsSec = startsSec(ichunk);
stopsSec = stopsSec(ichunk);

NRF = length(startsSec);
disp(['Number of data chunks to evaluate ',num2str(NRF)])

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
    disp(['Detecting in chunk ',num2str(k),' in ',fname])
    %     fid = fopen(fullfile(dname,fname(R,:)),'r');
    fid = fopen(fname,'r');
    ts = [];
    ts = ioReadXwav(fid,hdr,startsSec(k),stopsSec(k),channel,ftype);
    %     skip = byte_loc(R);
    %     fseek(fid,skip,-1);    % skip to correct rawfile
    spts = hdr.xhd.byte_length(k)/bps;
    %     ts = fread(fid,spts,'int16');   % read one time slice ie one raw file
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
            T(c) = datenum([0 0 0 0 0 (tds + startsSec(k))]) + hdr.start.dnum; % max amplitude detection datenum time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 140919 smw add timeseries snippet and spectra calc output
            
            if J(c) <= dur/2 +2  % put zeros on the front (ie beginning of file)
                dc = dur/2 +2 - J(c);
                if dc > 0
                    snip(c,:) = [zeros(1,dc),fts(1:J(c)+dur/2)];
                    ufsnip(c,:) = [zeros(1,dc),ts(1:J(c)+dur/2)];
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
                    snip(c,:) = [fts(J(c)-dur/2-1:spts),zeros(1,dc)];
                    ufsnip(c,:) = [ts(J(c)-dur/2-1:spts),zeros(1,dc)];
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
    disp(['Number of Detections in file segment = ',num2str(Ndet)])
    %
    if (Ndet > 0)
        ppdb = 20.*log10(pp);     % logify
        %     snipdb = 20.*log10(snip);
        % apply transfer function
        specdb = 10.*log10(spec*fs/2);
        ufspecdb = 10.*log10(ufspec*fs/2);
        nspecdb = 10.*log10(nspec*fs/2);
        Ptfx = interp1(tfFreq,tfPower,f,'linear','extrap');
        %         for i = 1:Ndet
        specdb = specdb + ones(Ndet,1);
        ufspecdb = ufspecdb + ones(Ndet,1);
        nspecdb = nspecdb + ones(Ndet,1);
        %         end
        % group for output
        TT = [TT, T];
        PP = [PP, ppdb];
        SN = [SN; snip];
        USN = [USN; ufsnip];
        SP = [SP; specdb];
        USP = [USP; ufspecdb];
        NSP = [NSP; nspecdb];
    end
    fclose(fid);    % close raw file
end
disp(' ')

