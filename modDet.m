% modDetKO.m
% modify detections based on analyst input
% originally based on jah/smw GoM beaked whale detector work
% JAH
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
stn = input('Enter Site name (MC GC DT): ','s'); % site name
dpn = input('Enter Deployment number (01 02 ...): ','s');     % deployment number
itnum = input('Enter Iteration number (1 2 ...) IN: ','s');  
itnumo = input('Enter Iteration number (1 2 ...) OUT: ','s');  
sp = input('Enter Species: Zc Me BWG Md Ko De ','s');
if (strcmp(sp,'Ko') || strcmp(sp,'k'))
    specchar = 'K'; %Simone abbreviation for species
    spe = 'Kogia';  tfselect = 80000; % freq used for transfer function
    icimin = 0.05; icimax = .3;
    thres = 116; % dB threshold
elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
    specchar = 'Z'; %Simone abbreviations for BW species
    spe = 'Cuviers';  tfselect = 40200; % freq used for transfer function
    icimin = 0.05; icimax = .7;
    thres = 121; % dB threshold
elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
    specchar = 'M'; %Simone abbreviations for BW species
    spe = 'Gervais'; tfselect = 40200; % freq used for transfer function
    icimin = 0.05; icimax = .5;
    thres = 121; % dB threshold
elseif (strcmp(sp,'BWG') || strcmp(sp,'g'))
    specchar = 'G'; %Simone abbreviations for BW species
    spe = 'BWG';  tfselect = 40200; % freq used for transfer function
elseif (strcmp(sp,'Md') || strcmp(sp,'d'))
    specchar = 'D'; %Simone abbreviations for BW species
    spe = 'BW31';  tfselect = 40200; % freq used for transfer function
elseif (strcmp(sp,'De') || strcmp(sp,'de'))
    %specchar = 'D'; %Simone abbreviations for BW species
    spe = 'Delphin';  tfselect = 0; % freq used for transfer function
    icimin = 0.05; icimax = .3;
    thres = 136.9; % dB threshold
else
    disp(' Bad Species type')
    return
end
%  Transfer Function
disp('Load Transfer Function File');
[fname,pname] = uigetfile('I:\Harp_TF\*.tf','Load TF File');
tffn = fullfile(pname,fname);
if strcmp(num2str(fname),'0')
    disp('Cancelled TF File');
    return
else %give feedback
    disp(['TF File: ',tffn]);
end
fid = fopen(tffn);
[A,count] = fscanf(fid,'%f %f',[2,inf]);
tffreq = A(1,:);
tfuppc = A(2,:);
fclose(fid);
if (tfselect > 0)
    tf = interp1(tffreq,tfuppc,tfselect,'linear','extrap');
    disp(['TF @',num2str(tfselect),' Hz =',num2str(tf)]);
else
    tf = 0;
    disp('No TF Applied');
end
% detection file
disp('Select Directory with Detections');
sdir = uigetdir('I:\','Select Directory with Detections');
pn1 = [sdir,'\'];
dpst = [stn,dpn];
% Current Detections
MTT = [];  MPP = [];  MSN = [];  MSP =[];
ia = []; ic = [];
fn1 = [dpst,'_',spe,'_TPWS',itnum,'.mat'];
DTfn = fullfile(pn1,fn1);% detection time file
load(DTfn)  % detection time and RL vectors: MTT, MPP, MSN, MSP
[DT0,ia,ic] = unique(MTT);   
disp([' non-unique removed: ',num2str(length(ic)- length(ia))]);
RL0 = MPP(ia)' + tf;
SN0 = MSN(ia,:);  % n by 202
SP0 = MSP(ia,:);  % n by 101
disp([' DT0: ',num2str(length(DT0))]);
% remove low amplitude 
ib = find(RL0 >= thres);
disp([' Removed too low:',num2str(length(ia)-length(ib))]);
DTb = DT0(ib);
RLb = RL0(ib)';
SNb = SN0(ib,:);  
SPb = SP0(ib,:);  
% False Detections
fn2 = [dpst,'_',spe,'_FD',itnum,'.mat'];
FDfn = fullfile(pn1,fn2); % false detection time file
load(FDfn)  % false detections vector : zFD
FD = unique(zFD);
disp([' FD: ',num2str(length(FD)),' Current: ', ...
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
MTT = DT1'; MPP = RL1' - tf; MSN = SN1; MSP = SP1;
fn8 = [dpst,'_',spe,'_TPWS',itnumo,'.mat'];
ofn = fullfile(pn1,fn8);
save(ofn,'MTT','MPP','MSN','MSP');
disp('Done Modifying File')
% Make ICI and PP plots and 5 min bins
Calicippfunc(MTT,RL1',stn,dpn,spe,pn1,icimin,icimax)
%binDur bin duration = 1 minute
%CalBinfun(stn,dpn,itnumo,1)
%



