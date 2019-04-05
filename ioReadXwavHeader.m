function xhd = ioReadXwavHeader(fileName, DateRE)
% 3/2/2011
% rdxwavhd.m
%
% reads pseudo-wav (XWAV or *.x.wav) file header
%
% functionized it for triton (less general, but puts values in global
% varibable PARAMS space
% smw 20 Oct, 2004
%
% smw 3-12 Aug, 2004 update again...introduced gain
%
% 060203smw updated for using all timing headers
%
% 060610 smw renamed xhd.SubchunkID and SubchunkSize with
% prefixes f,h,d for format, harp, and data subchunks
%
%clear all
%clc

% make xhd empty
xhd = [];

fid =  fopen(fileName,'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RIFF chunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhd.ChunkID = char(fread(fid,4,'uchar'))';       % "RIFF"
xhd.ChunkSize = fread(fid,1,'uint32');           % File size - 8 bytes
fileSize = getfield(dir(fileName),'bytes');
if xhd.ChunkSize ~= fileSize - 8
    disp('Error - incorrect Chunk Size')  
    % return    % comment to work with bad files
end
xhd.Format = char(fread(fid,4,'uchar'))';        % "WAVE"

if ~strcmp(xhd.ChunkID,'RIFF') || ~strcmp(xhd.Format,'WAVE')
    disp('not wav file - exit')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format Subchunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhd.fSubchunkID = char(fread(fid,4,'uchar'))';    % "fmt "
xhd.fSubchunkSize = fread(fid,1,'uint32');        % (Size of Subchunk - 8) = 16 bytes (PCM)
xhd.AudioFormat = fread(fid,1,'uint16');         % Compression code (PCM = 1)
xhd.NumChannels = fread(fid,1,'uint16');         % Number of Channels
xhd.SampleRate = fread(fid,1,'uint32');          % Sampling Rate (samples/second)
xhd.ByteRate = fread(fid,1,'uint32');            % Byte Rate = SampleRate * NumChannels * BitsPerSample / 8
xhd.BlockAlign = fread(fid,1,'uint16');          % # of Bytes per Sample Slice = NumChannels * BitsPerSample / 8
xhd.BitsPerSample = fread(fid,1,'uint16');       % # of Bits per Sample : 8bit = 8, 16bit = 16, etc

if ~strcmp(xhd.fSubchunkID,'fmt ') || xhd.fSubchunkSize ~= 16
    disp_msg('unknown wav format - exit')
    return
end

% should only be needed for special case bad data
% remove after debugged
% if xhd.SampleRate == 100000
%     disp_msg('Warning, changing sample rate from 100,000 to 500,000 Hz')
% %     xhd.SampleRate = 200000;
% %     xhd.ByteRate = 400000;
%     xhd.SampleRate = 500000;
%     xhd.ByteRate = 1000000;
% end

% copy to another name, and get number of bytes per sample
nBits = xhd.BitsPerSample;       % # of Bits per Sample : 8bit = 8, 16bit = 16, etc
% samp.byte = floor(nBits/8);       % # of Bytes per Sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARP Subchunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhd.hSubchunkID = char(fread(fid,4,'uchar'))';    % "harp"
if strcmp(xhd.hSubchunkID,'data')
    fs = xhd.SampleRate;
    disp('normal wav file - read data now')
    return
elseif ~strcmp(xhd.hSubchunkID,'harp')
    disp('unsupported wav format')
    disp(['SubchunkID = ',xhd.hSubchunkID])
    return
end
xhd.hSubchunkSize = fread(fid,1,'uint32');        % (Size of Subchunk - 8) includes write subchunk
xhd.WavVersionNumber = fread(fid,1,'uchar');     % Version number of the "harp" header (0-255)
xhd.FirmwareVersionNumber = char(fread(fid,10,'uchar'))';  % HARP Firmware Vesion
xhd.InstrumentID = char(fread(fid,4,'uchar'))';         % Instrument ID Number (0-255)
xhd.SiteName = char(fread(fid,4,'uchar'))';             % Site Name, 4 alpha-numeric characters
xhd.ExperimentName = char(fread(fid,8,'uchar'))';       % Experiment Name
xhd.DiskSequenceNumber = fread(fid,1,'uchar');   % Disk Sequence Number (1-16)
xhd.DiskSerialNumber = char(fread(fid,8,'uchar'))';     % Disk Serial Number
xhd.NumOfRawFiles = fread(fid,1,'uint16');         % Number of RawFiles in XWAV file
xhd.Longitude = fread(fid,1,'int32');           % Longitude (+/- 180 degrees) * 100,000
xhd.Latitude = fread(fid,1,'int32');            % Latitude (+/- 90 degrees) * 100,000
xhd.Depth = fread(fid,1,'int16');               % Depth, positive == down
if xhd.WavVersionNumber == 2
    xhd.drate = fread(fid,xhd.NumChannels,'single');
end
xhd.Reserved = fread(fid,8,'uchar')';            % Padding to extend subchunk to 64 bytes

if xhd.WavVersionNumber == 2
    hscs = (64 + 4*xhd.NumChannels) - 8 + xhd.NumOfRawFiles * (32 + 4*xhd.NumChannels);
else
    hscs = 64 - 8 + xhd.NumOfRawFiles * 32;
end
if xhd.hSubchunkSize ~= hscs
    sprintf('Error - HARP SubchunkSize and NumOfRawFiles discrepancy?\n')
    sprintf(['hSubchunkSize = %s\n',num2str(xhd.hSubchunkSize)])
    sprintf(['NumOfRawFiles = %s\n',num2str(xhd.NumOfRawFiles)])
 %   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write sub-sub chunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:xhd.NumOfRawFiles
                                                        % Start of Raw file :
    xhd.year(i) = fread(fid,1,'uchar');          % Year
    xhd.month(i) = fread(fid,1,'uchar');         % Month
    xhd.day(i) = fread(fid,1,'uchar');           % Day
    xhd.hour(i) = fread(fid,1,'uchar');          % Hour
    xhd.minute(i) = fread(fid,1,'uchar');        % Minute
    xhd.secs(i) = fread(fid,1,'uchar');          % Seconds
    xhd.ticks(i) = fread(fid,1,'uint16');        % Milliseconds
    xhd.byte_loc(i) = fread(fid,1,'uint32');     % Byte location in xwav file of RawFile start
    xhd.byte_length(i) = fread(fid,1,'uint32');    % Byte length of RawFile in xwav file
    xhd.write_length(i) = fread(fid,1,'uint32'); % # of blocks in RawFile length (default = 60000)
    xhd.sample_rate(i) = fread(fid,1,'uint32');  % sample rate of this RawFile
    xhd.gain(i) = fread(fid,1,'uint8');          % gain (1 = no change)
    xhd.padding = fread(fid,7,'uchar');    % Padding to make it 32 bytes...misc info can be added here
    if xhd.WavVersionNumber == 2
        xhd.dt = fread(fid,xhd.NumChannels,'single');     % time diff to channel 1
    end
    
    % should only be needed for special case bad data
    % remove after debugging
%     if xhd.sample_rate(i) == 100000
%         %   disp('Warning, changing sample rate from 100,000 to 200,000 Hz')
% %         xhd.sample_rate(i) = 200000;
%          xhd.sample_rate(i) = 500000;
% 
%     end
    
%     % calculate starting time [dnum => datenum in days] for each raw
%     % write/buffer flush
%     raw.dnumStart(i) = datenum([xhd.year(i) xhd.month(i)...
%         xhd.day(i) xhd.hour(i) xhd.minute(i) ...
%         xhd.secs(i)+(xhd.ticks(i)/1000)]);
%     raw.dvecStart(i,:) = [xhd.year(i) xhd.month(i)...
%         xhd.day(i) xhd.hour(i) xhd.minute(i) ...
%         xhd.secs(i)+(xhd.ticks(i)/1000)];
%     
%     % end of RawFile:
%     raw.dnumEnd(i) = raw.dnumStart(i) ...
%         + datenum([0 0 0 0 0 (xhd.byte_length(i) - 2)  ./  xhd.ByteRate]);
%     raw.dvecEnd(i,:) = raw.dvecStart(i,:) ...
%         + [0 0 0 0 0 (xhd.byte_length(i) - 2)  ./  xhd.ByteRate];
%     raw.dnumEnd(i) = raw.dnumStart(i) ...
%         + datenum([0 0 0 0 0 ceil((xhd.byte_length(i))  ./  (xhd.ByteRate * xhd.NumChannels))]);
%     raw.dvecEnd(i,:) = raw.dvecStart(i,:) ...
%         + [0 0 0 0 0 (xhd.byte_length(i))  ./  (xhd.ByteRate * xhd.NumChannels)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA Subchunk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xhd.dSubchunkID = char(fread(fid,4,'uchar'))';    % "data"
if ~strcmp(xhd.dSubchunkID,'data')
    disp('hummm, should be "data" here?')
    fprintf(['SubchunkID = %0.0d\n',xhd.dSubchunkID])
    return
end
xhd.dSubchunkSize = fread(fid,1,'uint32');        % (Size of Subchunk - 8) includes write subchunk

% read some data and check
%data = fread(fid,[4,100],'int16');

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rename a few things
%
%nch = xhd.NumChannels;         % Number of Channels
%fs = xhd.sample_rate(1);      % Real xwav Sampling Rate (samples/second)
% fs = xhd.SampleRate;        % Sampling Rate(samples/second)
                                            % this is 'wav sample rate,
                                            % could be fake...

% vectors (NumOfWrites)
%xgain = xhd.gain;          % gain (1 = no change)
%samp = xhd.byte_length ./ (nch * samp.byte); % # of Samples = Bytes of data / (# of Bytes per Sample Slice = NumChannels * BitsPerSample / 8)

%start.dnum = raw.dnumStart(1);
%start.dvec = raw.dvecStart(1,:);
%.dnum = raw.dnumEnd(xhd.NumOfRawFiles);
end