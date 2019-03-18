function params = getParams(userParams,varargin)
%getParams get default and user defined parameters
%
% Get parameters for the interface
%
% Copyright(C) 2019 by John A. Hildebrand, UCSD, jahildebrand@ucsd.edu
%                      Kait E. Frasier, UCSD, krasier@ucsd.edu
%                      Alba Solsona Berga, UCSD, asolsonaberga@ucsd.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
spParams = initSpParams('sp',sp,'srate',srate);

% General parameters (directories, iterations)
params.filePrefix = filePrefix;
params.itnum = itnum;
params.srate = srate;
params.sdir = sdir;
params.tfName = tfName;

%% create struct to return parameters
switch analysis
    case {'detEdit','mkLTSA'} 
        params.speName = speName;
        params.tfSelect = tfSelect;
        params.dtHi = dtHi;
        params.fLow = fLow;
        params.fHi = fHi;
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
        params.dfManual = dfManual;
        params.dfManual = dfManual;
        params.minBout = minBout;
        params.minDur = minDur;
        params.slope = slope;
        params.binDur = binDur;
        params.rawFileDur = rawFileDur;
        params.c4fd = c4fd;
        params.specploton = specploton;
        params.minNdet = minNdet; 
        params.maxDetLoad = maxDetLoad;
        params.gth = gth;
        params.colorTab = colorTab;
        
        % apply species default parameters
        if exist('spParams','var')
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
        params.tfSelect = tfSelect;
        params.threshRL = threshRL;
        params.dbRange = dbRange;
        params.iciRange = iciRange;
        params.frRange = frRange;
        params.frdbwRange = frdbwRange;
        params.durRange = durRange;
        params.durstep = durstep;
        params.N = N;
        params.gth = gth;
        params.minBout = minBout;
        params.p1Hi = p1Hi;
        
        % apply species default parameters
        if exist('spParams','var')
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
        
    case {'SumPPICIBin'}
        params.iciRange = iciRange;
        params.speName = speName;
        params.threshRL = threshRL;
        params.gth = gth;
        params.p1Hi = p1Hi;
        params.binDur = binDur;
        
        % apply species default parameters
        if exist('spParams','var')
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

