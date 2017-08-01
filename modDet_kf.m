% modDetKO.m
% modify detections based on analyst input
% originally based on jah/smw GoM beaked whale detector work
%
% 140305 smw
% 140318 jah  141003 jah modified for TPWS files
% 150813 jah for KO
%
% algorithm:
% per site and species
% DT1 = DT0 - FD
% where
% DT0 = current true detection times
% DT1 = modifed detection times
% from analyst input :
% FD = false detections times (eg ships, others species)
%
clear all
% thres = 116;  % detection threshold 116 dB
% get user input and set up file names
inDir = 'I:\MC\MC01_TPWS';
outDir = 'I:\MC\MC01_TPWS2';
if ~isdir(outDir)
    mkdir(outDir)
end
fList = dir(fullfile(inDir,'*TPWS1.mat'));
itnumo = 2;
cd(inDir);
iciMin = 0.02; iciMax = .35;
thres = 120; % dB threshold
for iFile = 1:length(fList)
    inFile = fList(iFile).name;
    load(inFile)
    
    % detection file
    [DT0,ia,ic] = unique(MTT);
    disp([' non-unique removed: ',num2str(length(ic)- length(ia))]);
    RL0 = MPP(ia);
    SN0 = MSN(ia,:);  % n by 202
    SP0 = MSP(ia,:);  % n by 101
    disp([' DT0: ',num2str(length(DT0))]);
    % remove low amplitude
    ib = find(RL0 >= thres);
    disp([' Removed too low:',num2str(length(ia)-length(ib))]);
    DTb = DT0(ib);
    RLb = RL0(ib);
    SNb = SN0(ib,:);
    SPb = SP0(ib,:);
    % False Detections
    load(strrep(inFile,'TPWS1','FD1')) % false detections vector : zFD
    FD = unique(zFD);
    disp([' Previous FD: ',num2str(length(FD)),' Current FD: ', ...
        num2str(length(intersect(DTb,FD)))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    DTa = [];
    RLa = [];
    SNa = [];
    SPa = [];
    % remove false detections
    IA = [];
    [DTa,IA] = setdiff(DTb',FD); % data in DT0 that is not in FD
    RLa = RLb(IA);
    SNa = SNb(IA,:);
    SPa = SPb(IA,:);
    disp([' DTb - FD: ',num2str(length(DTa))]);
    % make unique and sort by time
    IC = []; IA = [];
    DT1 = []; RL1 = []; SN1 = []; SP1 = [];
    [DT1,IA,IC] = unique(DTa,'sorted');
    nonu = length(IC)- length(IA);
    disp([' non-unique removed: ',num2str(nonu)]);
    RL1 = RLa(IA);
    SN1 = SNa(IA,:);
    SP1 = SPa(IA,:);
    disp([' Final Detections: ',num2str(length(DT1))]);
    % save modified detection file out
    MTT = []; MPP = []; MSN = [];  MSP = [];
    MTT = DT1'; MPP = RL1; MSN = SN1; MSP = SP1;
    outFileTPWS = strrep(inFile,'TPWS1','TPWS2');
    outFileName = fullfile(outDir,outFileTPWS);
    if exist('f','var')
        save(outFileName,'f','MTT','MPP','MSN','MSP','-v7.3')
    else
        warning('no frequency vector available')
        save(outFileName,'MTT','MPP','MSN','MSP','-v7.3')
    end
    disp('Done Modifying File')
    
end
% Make ICI and PP plots and 5 min bins
%Calicippfunc(MTT,RL1',stn,dpn,spe,pn1,icimin,icimax)
%binDur bin duration = 1 minute
%CalBinfun(stn,dpn,itnumo,1)
%



