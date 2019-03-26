%% Density_ClickCount

clearvars
close all

%% Parameters defined by user
filePrefix = {'GC'}; % File name to match. 
% File prefix should include site. 
sp = 'Pm'; % your species code
itnum = '2'; % which iteration you are looking for
tpwsPath = 'H:\newTPWS'; %directory of TPWS files
binsize = 5; % minutes per bin
refTime = '2010-04-01'; %reference time format 'yyyy-MM-dd'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% define subfolder that fit specified iteration
if itnum > 1
   for id = 2: str2num(itnum) % iternate id times according to itnum
       subfolder = ['TPWS',num2str(id)];
       tpwsPath = (fullfile(tpwsPath,subfolder));
   end
end

outDir = fullfile(tpwsPath,'Densities');
if ~isdir(outDir)
    disp(['Make new folder: ',outDir])
    mkdir(outDir)
end

myData = [];
for isite = 1:length(filePrefix)
%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:
p = sp_density_params('sp',sp,'site',filePrefix{isite});

FileParams = [filePrefix{isite},'_',sp,'_summaryData_forDensity.mat'];
load(fullfile(tpwsPath,FileParams))

% convert count of clicks to clicks per 5min bin
clickPerBin = weekTable.Count_Click;
% Click Method
Density = 1000*((clickPerBin.*(1-p.fpRate))./...
        (pi*(p.maxRadius_km^2)*p.pDet*(weekTable.Effort_Sec)*p.clickRate));
meanDensity = nanmean(Density);
coeffVar = repmat((p.fpRateCV^2 + p.pDetCV^2 + p.clickRateCV^2),length(Density),1) ;
varDensity = Density.^2.*coeffVar;
stdevDensity = sqrt(varDensity);

% MC
% if isite == 1
%     myData = [myData; datenum(weekTable.tbin(1:190)), Density(1:190),stdevDensity(1:190)];
% elseif isite == 2
%     myData = [myData; datenum(weekTable.tbin(191:end)), Density(191:end),stdevDensity(191:end)];
% end

% The others
myData = [myData; datenum(weekTable.tbin), Density,stdevDensity];
end

addData = [(datenum(refTime):7:myData(1,1))',nan(length(datenum(refTime):7:myData(1,1)),2)];
addDataEnd = [(myData(end,1):7:datenum('2017-12-31'))', nan(length(myData(end,1):7:datenum('2017-12-31')),2)];
myData = [addData;myData;addDataEnd];
outputFileName = ['GOM_',sp];
Deseason_data(myData,outDir,outputFileName,[],0,filePrefix{isite})

