function mkLtsa
%mkLtsa make long-term spectral averages (ltsa) from WAV/XWAV files in a 
%directory
%
% Syntax:
%   mkLtsa
%
% Make a LTSA file from selected audio files by reading audio headers. 
% A LTSA (Long-Term Spectral Average) is a compressed format to view
% large data sets. It is only compatible with WAV and XWAV files and 
% supports 4 filename formats:
%   1.  yymmdd-HHMMSS
%   2.  yymmdd_HHMMSS
%   3.  yyyymmdd_HHMMSS 
%   4.  yymmddHHMMSS
%   5.  yyyymmddTHHMMS
%
% Input prompt:
%   file format     File format(1 = WAVE, 2 = XWAV)
%   file directory  Path to the files to be processed
%   save direcotry  Path to save ltsa file
%   parameters      Parameters for generating ltsa files, time average
%                   length and frequency bin size
%
% Output:
%   ltsa file       
%
% Do not modify the following line, maintained by CVS
% Copyright(C)
% $Id: mk_ltsa.m,v 1.8 2008/11/15 17:08:06 mroch Exp $
%
% Modified 2018 by John A. Hildebrand, UCSD, jahildebrand@ucsd.edu
%                  Kait E. Frasier, UCSD, krasier@ucsd.edu
%                  Alba Solsona Berga, UCSD, asolsonaberga@ucsd.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global PARAMS

% initialize ltsa parameters to seet the defaults first
initLtsaParams

% get directory of wave/xwav files
% prompt user for file format, (1 = WAVE, 2 = XWAV) 
ioGetLtsaDir

if PARAMS.ltsa.gen == 0
    disp('Canceled making ltsa')
    return
end

% read data file headers
ioGetHeaders

if PARAMS.ltsa.gen == 0 % could not read wav file metadata 
    return
end

% get ltsa parameters from user (time avg. length and freq. bin size)
ioGetLtsaParams

% check some ltsa parameter and adjusts/gives suggestions of better params
ckLtsaParams

% setup lsta file header + directory listing
ioWriteLtsaHead

if PARAMS.ltsa.gen == 0
    disp('Canceled making ltsa')
    return
end
% calculated averaged spectra
calcLtsa

clear global