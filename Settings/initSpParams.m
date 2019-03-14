function spParams = initSpParams(varargin)
% Establish basic parameter settings, then update for species
% specific defaults, and user preferences.
% Pulled into subroutine kf 10/4/2016

% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'sp'
            sp = varargin{n+1}; n=n+2;
        case 'srate'
            srate = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

% Set parameters according to sp
%%% Delphinids
if  strcmpi(sp,'De')
    spParams.speName = 'Delphin';  
  
    % Bout parameters 
    spParams.threshRL = 118;
    spParams.c4fd = 1000; 
    % Panel LTSA and time series
    spParams.rlLow = spParams.threshRL - 6.9;
    spParams.rlHi = 190;
    spParams.dtHi = .6;
    % Panel Frequency spectra
    spParams.fLow = 10;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sperm whales

elseif (strcmp(sp,'PM') || strcmp(sp,'pm') || strcmp(sp,'Pm'))
    spParams.speName = 'Pm'; 
    
    % Bout parameters 
    spParams.minBout = 75;
    spParams.dfManual = 100;
    spParams.threshRL = 125;
    spParams.c4fd = 3000;
    % Panel LTSA and time series
    spParams.ltsaContrast = 180; 
    spParams.ltsaBright = 73;
    spParams.dtHi = 2;
    spParams.minDur = 60;
    % Panel Frequency spectra
    spParams.fLow = 5;
    % Panel RL rms vs. RL pp | Peak freq.
    spParams.slope = .7;
    spParams.threshRMS = 95;
    spParams.threshPP = 140;
    spParams.threshHiFreq = 30;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [300, 2000];
    spParams.N = 512;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beaked whales

elseif (strcmp(sp,'Zc') || strcmp(sp,'z') || strcmp(sp,'Cuviers'))
    spParams.speName = 'Cuviers'; 
    
    % Bout parameters 
    spParams.tfSelect = 40200; 
    spParams.dfManual = 100;
    spParams.threshRL = 121;
    spParams.c4fd = 3000;
    % Panel LTSA and time series
    spParams.ltsaContrast = 150; 
    spParams.ltsaBright = 50;
    spParams.dtHi = 1;
    spParams.minDur = 60;
    % Panel Frequency spectra
    spParams.fLow = 10;
    % Panel RL rms vs. RL pp | Peak freq.
    spParams.threshRMS = 55;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [40, 750];
    spParams.dbRange = [90, 170];
    spParams.durRange = [30, 300];
    spParams.frdbwRange = [0, 80];
    spParams.durstep = 2;  
    
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    spParams.speName = 'Gervais'; 
    
    % Bout parameters 
    spParams.tfSelect = 40200; 
    spParams.threshRL = 121;
    spParams.c4fd = 3000;
    % Panel LTSA and time series
    spParams.dtHi = 1;
    % Panel Frequency spectra
    spParams.fLow = 25;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [40, 750];
    spParams.dbRange = [90, 170];
    spParams.durRange = [30, 300];
    spParams.frdbwRange = [0, 80];
    spParams.durstep = 2; 
    
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    spParams.speName = 'BWG';
    
    % Bout parameters 
    spParams.tfSelect = 40200;
    spParams.threshRL = 121;
    spParams.c4fd = 3000;
    % Panel LTSA and time series
    spParams.dtHi = 1;
    % Panel Frequency spectra
    spParams.fLow = 25;
    
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    spParams.speName = 'BW31';
    
    % Bout parameters 
    spParams.tfSelect = 40200; 
    spParams.dfManual = 100;
    spParams.threshRL = 121; 
    spParams.c4fd = 3000;
    % Panel LTSA and time series
    spParams.ltsaContrast = 150; 
    spParams.ltsaBright = 50;
    spParams.dtHi = 1;
    spParams.minDur = 60;
    % Panel Frequency spectra
    spParams.fLow = 10;
    % Panel RL rms vs. RL pp | Peak freq.
    spParams.threshRMS = 55;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [40, 750];
    spParams.dbRange = [90, 170];
    spParams.durRange = [30, 300];
    spParams.frdbwRange = [0, 80];
    spParams.durstep = 2;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kogia

elseif (strcmp(sp,'Ko') || strcmp(sp,'k'))
    spParams.speName = 'Kogia';
    
    % Bout parameters 
    spParams.tfSelect = 80000;
    spParams.threshRL = 116;
    spParams.c4fd = 3000;
    % Panel Frequency spectra
    spParams.fLow = 70;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [40, 130];
    spParams.dbRange = [100, 150];
    spParams.durRange = [30, 111];
    spParams.frRange = [80, 160];
    spParams.frdbwRange = [0, 80];
    spParams.durstep = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Porpoise

elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
    spParams.speName = 'Porpoise'; 

    % Bout parameters 
    spParams.threshRL = 100;
    % Panel LTSA and time series
    spParams.rlLow = spParams.threshRL - 5;
    spParams.rlHi = 190;
    spParams.ltsaContrast = 250; 
    spParams.ltsaBright = 100;
    spParams.dtHi = .5;
    % Panel Frequency spectra
    spParams.fLow = 25;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Beluga

elseif strcmpi(sp,'Dl') 
    spParams.speName = 'Beluga'; 
    
    % Bout parameters 
    spParams.tfSelect = 45000;
    spParams.threshRL = 110;
    % Panel LTSA and time series
    spParams.rlLow = spParams.threshRL - 5;
    spParams.rlHi = 170;
    spParams.ltsaContrast = 200; 
    spParams.ltsaBright = 70;
    spParams.dtHi = .5;    
    % Panel Frequency spectra
    spParams.fLow = 20;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [20, 500];
    spParams.dbRange = [90, 170];
    spParams.durRange = [10, 300];
    spParams.durstep = 2; 
    
%%% Narwhal    
elseif (strcmp(sp,'Mm') || strcmp(sp,'Narwhal'))
    spParams.speName = 'Narwhal';

    % Bout parameters 
    spParams.tfSelect = 48000;
    spParams.threshRL = 100;
    % Panel LTSA and time series
    spParams.ltsaContrast = 240; 
    spParams.ltsaBright = 35;
    spParams.dtHi = 1;    
    % Panel Frequency spectra
    spParams.fLow = 2;
    % parameters for modDet plots (if specified to compute)
    spParams.iciRange = [40, 1000];
    spParams.durRange = [10, 300];
    spParams.rawFileDur = 300; 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Others

elseif strcmpi(sp,'MFA')
    spParams.speName = 'MFA';  

    % Bout parameters 
    spParams.tfSelect = 4000;
    spParams.dfManual = 10;
    spParams.threshRL = 80;
    % Panel LTSA and time series
    spParams.rlLow = spParams.threshRL - 5;
    spParams.rlHi = 180;
    spParams.ltsaMax = .5;
    spParams.dtHi = 2;
    % Panel Frequency spectra
    spParams.fLow = 2;
       
elseif strcmpi(sp,'whs')
    spParams.speName = 'whs';    

    % Bout parameters 
    spParams.dfManual = 10;
    % Panel LTSA and time series
    spParams.rlLow = 0;
    spParams.rlHi = 20;
    spParams.ltsaContrast = 310; 
    spParams.ltsaBright = 100;
    spParams.ltsaLims = [5,30];
    spParams.ltsaMax = .5;
    spParams.dtHi = 2;
    % Panel Frequency spectra
    spParams.fLow = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else
    warning('Unknown Species Type!!!')
end



