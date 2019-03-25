% make_TPWS.m

% Example script for converting a folder of mat files containing detections 
% into a single TPWS file using make_TPWS_vars.m

% Define input/output locations
inDir = 'E:\MyFolder'; % identify folder containing detection files.
myFileFlag = '*.mat'; % include a string to match for identifying the files to be processed.
iterationNum = 1; % Iteration number, usually 1 when creating TPWS files from scratch.
saveDir = 'E:\MyTPWSfolder';
outputFileName = 'SiteA_2012';

% Identify variables containing detection parameters 
% (See comments in make_TPWS_vars.m for more detail on expected inputs)

myTimeSeries = timeSeriesMat; % Replace "timeSeriesMat" with the name of the 
% variable containing detection time series in the detection files to be
% processed. REQUIRED.

myDetectionTimes = detTimesMat; % Replace "detTimesMat" with the name of the 
% variable containing a vector of detection times in the detection files to be
% processed. REQUIRED.

mySpectra = spectraMat; % Replace "spectraMat" with the name of the 
% variable containing detection spectra in the detection files to be
% processed. Optional, can be []. If NOT provided, 'fftLength' and 'sampleRate'
% are required.

f = freqVal; % Replace "freqVal" with the name of the variable containing 
% frequency vector associated with mySpectra. Required if mySpectra is provided.

myPeak2PeakAmp = ppVals; % Replace "ppVals" with the name of the 
% variable containing peak-to-peak recieved levels in the detection files to be
% processed. Optional. If not provided, myTfFreq and myTfFunValues is
% required to estimate this.

myFFTLength = 400; % Replace with your desired FFT or a variable name to 
% be loaded. Required if spectraMat is NOT provided.

mySampleRate = 200000; % Replace with your sample rate or a variable name
% be loaded. Required if spectraMat is NOT provided.

myTfFunValues =tfValuesVector; % Replace with the name of the variable 
% containing your transfer function values. Required if spectra must be
% calculated peak to recieved levels must be calculated. 

myTfFrequency = tfFreqVector; % Replace with the name of the variable 
% containing the frequencies associated with your transfer function values.
% Required if spectra must be calculated peak to recieved levels must be calculated. 

myBandpassEdges = [5,95];% replace with your bandpass filter values or 
% a variable. Optional, but highly recommended if spectra are being
% calculated. Removes roll off regions to improve spectral display.


%% Begin calculations
% Find all the files to be processed
detectionFileList = dir(fullfile(inDir,myFileFlag));

% Initialize variables
fnew = [];
MSP = [];
MSN = [];
MTT = [];
MPP = [];

% Iterate over each file, collecting and calculating the various
% parameters.
for iFile = 1:length(detectionFileList)
    % load each detection file
    thisDetFile = fullfile(detectionFileList(iFile).folder,...
        detectionFileList(iFile).name);
    load(thisDetFile)
    
    [MPPnew,MTTnew,MSPnew,MSNnew,fnew] = make_TPWS_vars('timeSeries',myTimeSeries,...
        'detectionTimes',myDetectionTimes,'signalSpectra',mySpectra,'f',f,...
        'fftLength',myFFTLength,'peak2peakAmp',myPeak2PeakAmp,...
        'sampleRate',mySampleRate,'tfFunFrequency',myTfFrequency,...
        'tfFunValues',myTfFunValues,'bandPassEdges',myBandpassEdges);
    
    MPP = [MPP;MPPnew];
    MTT = [MTT;MTTnew];
    MSN = [MSN;MSNnew];
    MSP = [MSP;MSPnew];
    
    if ~isempty(f)
        f = fnew;
    end
    % Optionally can add logic in here for rolling into a new output file
    % if things get too large (>1.5 million detections is a rough estimate 
    % of max desirable file size). Preallocating space for MSP and MSN will
    % improve performance.
end

% Save output
fullOutputFileName = [outputFileName,'_TPWS',num2str(iterationNum),'.mat'];
save(fullfile(saveDir,fullOutputFileName),'MPP','MSP','MSN','MTT','f','-v7.3')

