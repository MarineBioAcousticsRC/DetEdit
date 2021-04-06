function params = getParams(userParams,varargin)

% getParams 
%
% Read initial and species default parameters, overwrite with user-defined 
%
% Inputs:
%
%   userFunc - Script user parameter settings.
%
%   'analysis' - REQUIRED.
%       A string specifying the type of analysis to get the specific
%       parameters.
%       Three options are accepted:
%           'detEdit' - interface parameters
%           'mkLTSA' - parameters to calculate the LTSA sessions
%           'modDet' - parameters to modify annotation files
%           'SumPPICIBin' - parameters to summarize  
%
%
% Output:
%
%   params - A struct with variable fields of parameter settings



% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'analysis'
            analysis = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

% Initialize default parameters first
initDefaultParams

% User defined parameters
userParams()

% Initialize defined species parameters
spParams = initSpParams('sp',sp,'sampleRate',sampleRate);

% General parameters (directories, iterations)
params.filePrefix = filePrefix;
params.iterationNum = iterationNum;
params.sampleRate = sampleRate;
params.tpwsDir = tpwsDir;
if exist('IDfilePrefix')
    if ~isempty(IDfilePrefix)
        params.IDfilePrefix = IDfilePrefix;
    else
        params.IDfilePrefix = filePrefix;
    end
end
params.tfName = tfName;
params.ltsaDir = ltsaDir;

%% create struct to return parameters
switch analysis
    case {'detEdit','mkLTSA'} 
        params.speName = speName;
        params.tfSelect = tfSelect;
        params.dtHi = dtHi;
        params.fLow = fLow;
        params.fHi = params.sampleRate/2;
        params.threshRL = threshRL;
        params.threshRMS = threshRMS;
        params.threshPP = threshPP;
        params.threshHiFreq = threshHiFreq;
        params.ltsaContrast = ltsaContrast;
        params.ltsaBright = ltsaBright;
        params.ltsaLims = ltsaLims;
        params.ltsaMax = ltsaMax;
        params.rlLow = rlLow;
        params.rlHi = rlHi;
        params.rmsLow = rmsLow;
        params.rmsHi = rmsHi;
        params.dfManual = dfManual;
        params.dfManual = dfManual;
        params.minBout = minBout;
        params.minDur = minDur;
        params.slope = slope;
        params.binDur = binDur;
        params.rawFileDur = rawFileDur;
        params.c4fd = c4fd; 
        params.nTestBins = [];
        params.specploton = specploton;
        params.minNdet = minNdet; 
        params.maxDetLoad = maxDetLoad;
        params.gth = gth;
        params.colorTab = colorTab;
        params.ltsaDir = ltsaDir;
        params.mySpID = mySpID;
        params.autoFalse = autoFalse;
        params.minLabelConfidence = minLabelConfidence;
        params.minClicks = minClicks;
        % apply species default parameters
        if ~isempty(spParams)
            idx = ismember(fieldnames(spParams),fieldnames(params));
            fn = fieldnames(spParams);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = spParams.(fnSel{f});
            end
        end
        
        % apply user parameters (it overwrites species default parameters)
        if exist('paramsUser','var')
            idx = ismember(fieldnames(paramsUser),fieldnames(params));
            fn = fieldnames(paramsUser);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = paramsUser.(fnSel{f});
            end
        end
        
    case {'modDet'} 
        params.speName = speName;
        params.tfSelect = tfSelect;
        params.threshRL = threshRL;
        params.dbRange = dbRange;
        params.iciRange = iciRange;
        params.frRange = frRange;
        params.frdbwRange = frdbwRange;
        params.durRange = durRange;
        params.durstep = durstep;
        params.gth = gth;
        params.minBout = minBout;
        params.excludeID = excludeID;
        params.calcParams = calcParams;
        params.ltsaMax = ltsaMax;
        
        % apply species default parameters
        if ~isempty(spParams)
            idx = ismember(fieldnames(spParams),fieldnames(params));
            fn = fieldnames(spParams);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = spParams.(fnSel{f});
            end
        end
        
        % apply user parameters (it overwrites species default parameters)
        if exist('paramsUser','var')
            idx = ismember(fieldnames(paramsUser),fieldnames(params));
            fn = fieldnames(paramsUser);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = paramsUser.(fnSel{f});
            end
        end
        
    case {'summaryParams'}
        params.tfSelect = tfSelect;
        params.iciRange = iciRange;
        params.speName = speName;
        params.threshRL = threshRL;
        params.gth = gth;
        params.rlHi = rlHi;
        params.binDur = binDur;
        params.effortTimes = effortTimes;
        params.referenceTime = referenceTime;
        
        % apply species default parameters
        if ~isempty(spParams)
            idx = ismember(fieldnames(spParams),fieldnames(params));
            fn = fieldnames(spParams);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = spParams.(fnSel{f});
            end
        end
        
        % apply user parameters (it overwrites species default parameters)
        if exist('paramsUser','var')
            idx = ismember(fieldnames(paramsUser),fieldnames(params));
            fn = fieldnames(paramsUser);
            fnSel = fn(idx);
            for f = 1:length(fnSel)
                params.(fnSel{f}) = paramsUser.(fnSel{f});
            end
        end
        
    otherwise
        sprintf(['No analysis specified. Please add one of these options:\n',...
        'detEdit, mkLTSA or modDet'])
end

