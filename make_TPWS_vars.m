function [MPP,MTT,MSP,MSN,f] = make_TPWS_vars(varargin)

% make_TPWS_vars.m

% Simple script takes lists of detection parameters and outputs TPWS
% variables. See make_TPWS for an example of combining output variables
% from multiple detection files into one TPWS file

% Inputs
%   'timeSeries' - REQUIRED. 
%       A 2D matrix where each row represents the time series of one detection. 
%       Alignment by maximum amplitude is recommended.
%
%   'detectionTimes' - REQUIRED
%       An Nx1 vector of detection times as Matlab datenums, where N is the
%       number of detection times.
%   
%   'signalSpectra' - Optional, will be calculated from time series if not
%       provided. 
%       An [N x nF] matrix where N is the number of detections and nF is
%       the number of frequency bins in the spectra.
%       If not provided, 'fftLength' and 'sampleRate' are required.
% 
%   'f' - Required if 'signalSpectra' is provided.
%       If not provided, will be calculated from 'fftLength' and 'sampleRate'.
%       A [1 x nF] vector of frequencies associated with 'spectra'. 
%       Units = kHz.
%
%   'fftLength' - Required if 'signalSpectra' is NOT provided. 
%       Number of points to use when calculating the FFT. 
% 
%   'sampleRate' - Required if 'signalSpectra' is NOT provided. 
%       Used to calculate frequency vector associated with spectra.
%       Units = Hz
% 
%   'peak2peakAmp' - Optional, will be calculated from time series if not
%       provided.
%       If NOT provided, a transfer function is required. 
%       Units: dB peak-to-peak, re: 1 \muPa.
%
%   'tfFunFrequency' - Required if 'signalSpectra' is NOT provided.
%                      Required if 'peak2peakAmp' is NOT provided, unless
%                      size('tfFunVals) = [1,1] (see 'tfFunValues' for details).
%       Units = kHz
%       A [1 x nTF] vector of frequencies associated with the transfer
%       function vector. Will be used to estimate peak to peak amplitude of
%       detections.
%
%   'tfFunValues' - Required if 'peak2peakAmp' is NOT provided.
%       Units = dB peak-to-peak, re: 1 \muPa.
%       Used to estimate peak to peak amplitude of detections
%       Two options are accepted:
%       a) A single number. In this case, the same adjustment is applied to 
%       all detections across all frequencies, and regardless of peak frequency.
%       b) A [1 x nTF] vector of amplitudes associated with the transfer
%       function frequency vector. Will be used to estimate peak-to-peak 
%       amplitude of detections, and adjust spectra according to sensor
%       sensitivity.
% 
%   'bandPassEdges' - Optional. Only used if signalSpectra is not provided
%       A [1x2] vector of upper and lower frequencies used to truncate the 
%       spectra if the input timeseries have been band passed. Values in
%       kHz.
%       
%       
%
% Output:
% 
%   MTT - An [N x 1] vector of detection times, where N is the
%         number of detections.
%   MPP - An [N x 1] vector of peak to peak received level (RL) amplitudes.
%   
%   MSP - An [N x nF] matrix of detection spectra, where nF is the length of
%        the spectra.
%   MSN - An [N x maxSampleWindow] matrix of detection timeseries.
%
%   f = An Fx1 frequency vector associated with MSPs

timeSeries = [];
detectionTimes = [];
signalSpectra = [];
f = [];
fftLength = [];
sampleRate = [];
peak2peakAmp = [];
tfFunFrequency = [];
tfFunValues = [];
bandPassEdges = [];

tfFunResampled = [];

vIdx = 1;
while vIdx <= length(varargin)
    switch varargin{vIdx}
        case 'timeSeries'
            timeSeries = varargin{vIdx+1};
            vIdx = vIdx+2;
        case 'detectionTimes'
            detectionTimes = varargin{vIdx+1};
            vIdx = vIdx+2;
        case 'signalSpectra'
            signalSpectra = varargin{vIdx+1};
            vIdx = vIdx+2;   
        case 'f'
            f = varargin{vIdx+1};
            vIdx = vIdx+2;   
        case 'fftLength'
            fftLength = varargin{vIdx+1};
            vIdx = vIdx+2;   
        case 'sampleRate'  
            sampleRate = varargin{vIdx+1};
            vIdx = vIdx+2;   
        case 'peak2peakAmp'
            peak2peakAmp = varargin{vIdx+1};
            vIdx = vIdx+2;      
        case 'tfFunFrequency'
            tfFunFrequency = varargin{vIdx+1};
            vIdx = vIdx+2;  
        case 'tfFunValues'
            tfFunValues = varargin{vIdx+1};
            vIdx = vIdx+2;  
        case 'bandPassEdges'
            bandPassEdges = varargin{vIdx+1};
            vIdx = vIdx+2; 
        otherwise
            sprintf('undefined input value %s',varargin{vIdx})
            vIdx = vIdx+2; 
    end
end
        
if isempty(timeSeries)
    error('timeSeries input is required')
elseif isempty(detectionTimes)
    error('detectionTimes input is required')
elseif isempty(f) && ~isempty(signalSpectra)
    error('f input is required because spectra were provided')
elseif isempty(fftLength) && isempty(signalSpectra)
    error('fftLength input is required to calculate spectra')
elseif isempty(sampleRate) && isempty(signalSpectra)
    error('sampleRate input is required to calculate spectra')
elseif isempty(tfFunValues) && (isempty(peak2peakAmp) || isempty(signalSpectra))
    error('tfFunValues input required to compute spectra & peak-to-peak amplitudes')
elseif isempty(tfFunFrequency) && length(tfFunValues) >1
    error('tfFunFrequency input required to go with tfFunValues vector')
end
% Make MSN
MSN = timeSeries;

% make MTT 
MTT = detectionTimes;

% Make MSP
if ~isempty(signalSpectra)
    MSP = signalSpectra;
else
    
    % Calculate window
    fftWindow = hann(fftLength)';
    
    % Pad waveforms with zeros to reach desired fft length
    padZeros = [MSN,zeros(size(MSN,1),fftLength-size(MSN,2))];
    
    % Compute power spectral density
    [MSP,fHz] = pwelch(padZeros',fftWindow,0,fftLength,sampleRate);

    % convert to dBs
    MSP = 10.*log10(abs(MSP'));

    f = fHz./1000; % Convert to kHz
    
    % Resample transfer function so it matches your frequency vector
    tfFunResampled = interp1(tfFunFrequency, tfFunValues, f, 'linear', 'extrap')';
    
    % Add transfer function to spectra
    MSP = MSP + repmat(tfFunResampled,size(MSP,1),1); 
    
    % Truncate spectra if filter edges are provided
    if ~isempty(bandPassEdges)
        [~,lowIdx] = min(abs(f-bandPassEdges(1)));
        [~,hiIdx] = min(abs(f-bandPassEdges(2)));
        MSP = MSP(:,lowIdx:hiIdx);
        f = f(lowIdx:hiIdx);
    end
end
    
% Compute MPP
if ~isempty(peak2peakAmp)
    MPP = peak2peakAmp;
else
    % Figure out the peak frequency of the spectra
    [~, peakFrIdx] = max(MSP,[],2);
    % Identify transfer fun value to add
    if ~isempty(tfFunResampled)
        % if we haven't done this above, do it here
    	tfFunResampled = interp1(tfFunFrequency, tfFunValues, f, 'linear', 'extrap')';
    end
    peakFrTfVal = tfFunResampled(peakFrIdx);
    % Compute amplitude from waveform and add tf.
    MPP = 20*log10(max(MSN,[],2)+abs(min(MSN,[],2))) + peakFrTfVal';
    
end

