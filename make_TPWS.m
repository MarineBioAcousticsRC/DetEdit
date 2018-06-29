% make_TPWS.m

% Example script for converting a folder of mat files containing detections 
% into a single TPWS file.

inDir = 'E:\MyFolder'; % identify folder containing detection files.
myFileFlag = '*.mat'; % include a string to match for identifying the files to be processed.
iterationNum = 1; %Iteration number, usually 1 when creating TPWS files from scratch.
saveDir = 'myTPWSfolder';

% Find all the files to be processed
detectionFileList = dir(fullfile(inDir,myFileFlag));

% Initialize variables
f = [];
MSP = [];
MSN = [];
MTT = [];
MPP = [];

% Iterate over each file, collecting and calculating the various
% parameters.
for iFile = 1:length(detectionFileList)
    thisDetFile = fullfile(detectionFileList(iFile).folder,...
        detectionFileList(iFile).name);
    load(thisDetFile)
    
    [MPPnew,MTTnew,MSPnew,MSNnew,fnew] = make_TPWS_vars('timeSeries',myTimeSeries,...
        'detectionTimes',myDetectionTimes,'fftLength',2000,'sampleRate',200000,...
        'tfFunFrequency', myTfFreq,'tfFunValues',myTfFunValues,...
        'bandPassEdges',[5,95]);
    
    MPP = [MPP;MPPnew];
    MTT = [MTT;MTTnew];
    MSN = [MSN;MSNnew];
    MSP = [MSP;MSPnew];
    
    if ~isempty(f)
        f = fnew;
    end
    % Optionally can put logic in here for rolling into a new output file
    % if things get too large.
end

% Save output
outputFileName = fullfile(saveDir,outputFileName);
save(fullfile(saveDir,outputFileName),'MPP','MSP','MSN','MTT','f')

