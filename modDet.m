function modDet(varargin)
% modDet.m
% modify detections based on analyst input
% originally based on jah/smw GoM beaked whale detector work
%
% 140305 smw
% 140318 jah  141003 jah modified for TPWS files
% 150813 jah for KO
%
% algorithm:
% per site and species
% DT1 = DT0 - zFD
% where
% DT0 = current true detection times
% DT1 = modifed detection times
% from analyst input :
% zFD = false detections times (eg ships, others species)

% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'filePrefix'
            filePrefix = varargin{n+1}; n=n+2;
        case 'detfn'
            detfn = varargin{n+1}; n=n+2;
        case 'sp'
            sp = varargin{n+1}; n=n+2;
        case 'sdir'
            sdir = varargin{n+1}; n=n+2;
        case 'srate'
            srate = varargin{n+1}; n=n+2;
        case 'itnum'
            itnum = varargin{n+1}; n=n+2;
        case 'getParams'
            getParams = varargin{n+1}; n=n+2;
        case 'tfName'
            tfName = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

secInDay = 60*60*24; % convert seconds to days

%% define output directory
newitnum = num2str(str2double(itnum)+1);
inTPWS = ['TPWS',itnum];
outTPWS = ['TPWS',newitnum];

outDir = fullfile(sdir,outTPWS);
if ~isdir(outDir)
    disp(['Make new folder: ',outDir])
    mkdir(outDir)
end

%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:
p = sp_setting_defaults('sp',sp,'srate',srate,'analysis','modDet');


%% user interface to get TF file
if (p.tfSelect > 0) || ~isempty(strfind(getParams,'all'))
    if ~exist('tfName','var')% user interface to get TF file
        disp('Load Transfer Function File');
        [fname,pname] = uigetfile('I:\Harp_TF\*.tf','Load TF File');
        tffn = fullfile(pname,fname);
    else % or get it automatically from tf directory provided in settings
        stndeploy = strsplit(filePrefix,'_'); % get only station and deployment
        tffn = findTfFile(tfName,stndeploy); % get corresponding tf file
    end
    disp(['TF File: ',tffn]);
    
    fid = fopen(tffn);
    [A,~] = fscanf(fid,'%f %f',[2,inf]);
    tffreq = A(1,:);
    tfuppc = A(2,:);
    fclose(fid);
    
    if(p.tfSelect > 0)
        tf = interp1(tffreq,tfuppc,p.tfSelect,'linear','extrap');
        disp(['TF @',num2str(p.tfSelect),' Hz =',num2str(tf)]);
    end
    if ~isempty(strfind(getParams,'all'))
        BinkHz = 0:1:srate/2;
        tf = interp1(tffreq,tfuppc,BinkHz,'linear','extrap');
        disp('TF Applied to get parameters');
    end
else
    tf = 0;
    disp('No TF Applied')
end

%% load original TPWS file
try
    load(fullfile(sdir,detfn))
catch
    error('Unable to read file %s. Not such file in directory: %s',detfn,sdir)
end

% inititalize
DT0 = []; DT1 = [];
RL0 = []; RL1 = [];
SN0 = []; SN1 = [];
SP0 = []; SP1 = [];

% define original detection files
DT0 = MTT;
RL0 = MPP;
SN0 = MSN;
SP0 = MSP;

% empty variables to fill later with modified detections
MTT = []; MPP = []; MSN = [];  MSP = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% False Detections
zFDfn = strrep(detfn,inTPWS,['FD',itnum]);
load(fullfile(sdir,zFDfn)) % false detections vector : zFD

% remove false detections
% DT1 = % data in DT0 that is not in zFD
[DT1,IA] = setdiff(DT0',zFD); % setdiff already sorts the data
RL1 = RL0(IA);
SN1 = SN0(IA,:);
SP1 = SP0(IA,:);
disp(['Number of Starting Detections = ',num2str(length(DT0))]);
disp(['Number of Final Detections = ',num2str(length(DT1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% ID Detections
zFDfn = strrep(detfn,inTPWS,['ID',itnum]);
load(fullfile(sdir,zFDfn)) % false detections vector : zFD

%%%%%%%%%%%%%%%%%%%%%%%%%%
% save modified detection to output file
outFileTPWS = strrep(detfn,inTPWS,outTPWS);
outFileID = strrep(outFileTPWS,'TPWS','ID');

MTT = DT1';
MSN = SN1;
MSP = SP1;
if (p.tfSelect > 0);
    MPP = RL1 + tf;
else
    MPP = RL1;
end

% check if there is at least one encounter longer than 75s, if not
% do not store TPWS file
if ~ isempty(MTT)
    dt = diff(MTT)*secInDay;
    I = find(dt>p.gth*60*60);
    sb = [MTT(1); MTT(I+1)];
    eb = [MTT(I); MTT(end)];
    bd = (eb - sb);
    bdI = find(bd > (p.minBout / secInDay)); % find bouts longer than the minimum (75s)
    
    if ~isempty(bdI)
        disp(['Save ',fullfile(outDir,outFileTPWS)])
        if exist('f','var')
            save(fullfile(outDir,outFileTPWS),'f','MTT','MPP','MSN','MSP','-v7.3')
        else
            warning('no frequency vector available')
            save(fullfile(outDir,outFileTPWS),'MTT','MPP','MSN','MSP','-v7.3')
        end
        disp(['Save ',fullfile(outDir,outFileID)])
        save(fullfile(outDir,outFileID),'zID','-v7.3')
        disp('Done Modifying File')
        
        % Calculate parameters and make figure if specified by user in itr_modDet
        switch getParams
            case 'ici&pp'
                Calicippfunc(MTT,MPP,filePrefix,sp,outDir,outFileTPWS,p)
            case 'all'
                CalPARAMSfunc(MTT,MPP',MSN,filePrefix,sp,outDir,outFileTPWS,p,tf,srate)
            case 'none'
                disp('No parameters calculated')
            otherwise
                fprintf(['Wrong name to call parameters, see itr_modDet:\n -all- for all parameters',...
                    '\n -ici&pp- for computing only ici and pp \n -none- No parameters'])
        end
    end
else
    disp(['No encounter longer than minimum bout (',num2str(p.minBout),') in file: ',outFileTPWS])
end
% add this in function Density
%binDur bin duration = 1 minute
%CalBinfun(stn,dpn,itnumo,1)




