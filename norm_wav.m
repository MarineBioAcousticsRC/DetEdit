function WAV = norm_wav(WAVSet)
% normalize waveform
mnWAV = min(WAVSet,[],2);
WAV = WAVSet - mnWAV;  % make low amp part = 0
mxWAV = max(WAV,[],2);
WAV = WAV ./ mxWAV;  % make high amp part = 1
