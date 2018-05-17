function initLtsaParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% init_ltsaparams.m
%
% initialize ltsa parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PARAMS

PARAMS.ltsa = [];
PARAMS.ltsahd = [];

% user defined:
PARAMS.ltsa.tave = 5;       % averaging time [seconds]
PARAMS.ltsa.dfreq = 200;    % frequency bin size [Hz]

% experiment defined (in XWAV file):
PARAMS.ltsa.fs = 200000;    % sample rate [Hz]

% calculated based on user defined:
% number of samples for fft (nfft = 1000 for dfreq=200Hz & fs=200000Hz)
PARAMS.ltsa.nfft = PARAMS.ltsa.fs / PARAMS.ltsa.dfreq;    
% compression factor (cfact = 1000 for tave=5sec,fs=200000Hz,dfreq=200)
PARAMS.ltsa.cfact = PARAMS.ltsa.tave * PARAMS.ltsa.fs / PARAMS.ltsa.nfft;   

% other
PARAMS.ltsa.indir = 'C:\';   % starting data directory
PARAMS.ltsa.ftype = 2;      % 1= WAVE, 2=XWAV

PARAMS.ltsa.dtype = 1;      % 1 = HARP, 2 = ARP, 3 = OBS, 4 = towed array or sonobuoy

PARAMS.ltsa.ch = 1;       % channel to do ltsa on for multichannel wav files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.ltsa.tseg.step = -1;     % Step size (== dur by default)
PARAMS.ltsa.tseg.hr = 2;
PARAMS.ltsa.tseg.sec = PARAMS.ltsa.tseg.hr * 60 * 60;         % initial window time segment duration

PARAMS.ltsa.ftype = 1;
PARAMS.ltsa.freq0 = 0;			% set frequency PARAMS.ltsa lower limit
PARAMS.ltsa.freq1 = -1;         % set frequency PARAMS.ltsa.ltsa upper limit
PARAMS.ltsa.bright = 0;			% shift in dB
PARAMS.ltsa.contrast = 100;		% amplify in % dB
PARAMS.ltsa.fax = 0;            % linear or log freq axis
PARAMS.ltsa.cmap = 'jet';		% color map for spectrogram
PARAMS.ltsa.start.yr = 0;
PARAMS.ltsa.start.str = '0000';
PARAMS.ltsa.aptime = 0;			%  pause time (typically CPU speed dependent?
PARAMS.ltsa.cancel = 0;
PARAMS.ltsa.delimit.value = 0;  %  delimit value is off at first


% PARAMS.ltsa.infile = PARAMS.ltsa.outfile;
% PARAMS.ltsa.inpath = PARAMS.ltsa.outdir;

