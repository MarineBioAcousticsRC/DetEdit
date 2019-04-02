function SumPPICIBin(userFunc)

%% Load Settings preferences
% Get parameter settings worked out between user preferences, defaults, and
% species-specific settings:

% get user input and set up function name
typeInput = exist('userFunc','var');
if typeInput ~= 1
    [userfile,userpath] = uigetfile('*.m',...
        'Select Script with your Data Parameter Settings');
    addpath(userpath) % it adds user folder path to the beggining of the set path
    userFunc = str2func(['@',userfile(1:end-2)]);
end

p = getParams(userFunc,'analysis','summaryParams');

%% Define subfolder that fit specified iteration
if p.iterationNum > 1
    for id = 2: str2num(p.iterationNum) % iternate id times according to p.iterationNum
        subfolder = ['TPWS',num2str(id)];
        p.tpwsDir = (fullfile(p.tpwsDir,subfolder));
    end
end

%% Check if TPWS file exists (does not look in subdirectories)
% Concatenate parts of file name
if isempty(p.speName)
    detfn = [p.filePrefix,'.*','TPWS',p.iterationNum,'.mat'];
else
    detfn = [p.filePrefix,'.*',p.speName,'.*TPWS',p.iterationNum,'.mat'];
end
% Get a list of all the files in the start directory
fileList = cellstr(ls(p.tpwsDir));
% Find the file name that matches the p.filePrefix
fileMatchIdx = find(~cellfun(@isempty,regexp(fileList,detfn))>0);
if isempty(fileMatchIdx)
    % if no matches, throw error
    error('No files matching filePrefix found!')
end

matchingFile = fileList{fileMatchIdx};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Handle Transfer Function
if (p.tfSelect > 0)
    [tf,~,~] = getTransfunc(p.filePrefix, p.tfName,p);
else
    tf = 0;
    disp('No TF Applied')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get effort times matching prefix file
allEfforts = readtable(p.effortTimes);
effTable = allEfforts(ismember(allEfforts.Sites,p.filePrefix),:);

% make Variable Names consistent
startVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'Start.*Effort'))>0,1,'first');
endVar = find(~cellfun(@isempty,regexp(effTable.Properties.VariableNames,'End.*Effort'))>0,1,'first');
effTable.Properties.VariableNames{startVar} = 'Start';
effTable.Properties.VariableNames{endVar} = 'End';

Start = datetime(x2mdate(effTable.Start),'ConvertFrom','datenum');
End = datetime(x2mdate(effTable.End),'ConvertFrom','datenum');

effort = timetable(Start,End);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Concatenate matrices
PPall = []; TTall = []; ICIall = []; % initialize matrices
for idsk = 1 : length(fileMatchIdx)

    % Find files and load data
    if iscell(matchingFile)
        detfn = dir(fullfile(p.tpwsDir,matchingFile{idsk}));
    else
        detfn = dir(fullfile(p.tpwsDir,matchingFile));
    end
    fprintf('Loading %d/%d file %s\n',idsk,length(fileMatchIdx),detfn.name)
    D = load(fullfile(p.tpwsDir,detfn.name));
    
    % find times outside effort (sometimes there are detections
    % which are from the audio test at the beggining of the wav file)
    within = cell2mat(arrayfun(@(x)sum(isbetween(x,datenum(effort.Start),datenum(effort.End))),D.MTT,'uni',false));
    goodIdx = find(within ~= 0);
    MTT = D.MTT(goodIdx); % only keep the good detections
    MPP = D.MPP(goodIdx) + tf;
    MSP = D.MSP(goodIdx,:);
    
    % concatenate
    TTall = [TTall; MTT];   % group start times
    PPall = [PPall; MPP];   % group peak-to-peak
    
    % Inter-Click Interval
    ici = diff(MTT)*24*60*60*1000; % in ms
    ICIall = [ICIall;[ici; nan]];  % group inter-click interval
end

%% After parfor data may not be sorted. Sort all the variables
[~,sorted] = sort(TTall);
TTall = TTall(sorted);
PPall = PPall(sorted);
ICIall = ICIall(sorted);

%% Create timetable per click
tbin = datetime(TTall,'ConvertFrom','datenum');
clickData = timetable(tbin,PPall,ICIall);
clear tbin
%% Convert times to bin vector times
vTT = datevec(TTall);
tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
    zeros(length(vTT),1)]);

%% create table and get click counts and max pp per bin
data = timetable(tbin,TTall,PPall);

binData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
binData.Properties.VariableNames{'GroupCount'} = 'Count'; % #clicks per bin
binData.Properties.VariableNames{'max_PPall'} = 'maxPP';

positiveCounts = sum(binData.Count);
positiveBins = length(binData.Count);

%% group effort in bins
effort.diffSec = seconds(effort.End-effort.Start) ;
effort.bins = effort.diffSec/(60*p.binDur);
effort.roundbin = round(effort.diffSec/(60*p.binDur));

% convert intervals in bins 
binEffort = intervalToBinTimetable(effort.Start,effort.End,p); 
binEffort.Properties.VariableNames{1} = 'bin';
binEffort.Properties.VariableNames{2} = 'sec';

secMonitEffort = sum(binEffort.sec);
binMonitEffort = sum(binEffort.bin);

%% get average of detection by effort
NktTkt = positiveCounts/secMonitEffort;
NktTktbin = positiveBins/binMonitEffort;
disp(['Nkt/Tkt for click counting: ',num2str(NktTkt)]);
disp(['Nkt/Tkt for group counting: ',num2str(NktTktbin)]);

%% group data by days and weeks
Click = retime(binData(:,1),'daily','sum'); % # mean click per day
Bin = retime(binData(:,1),'daily','count'); % #bin per day

dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';

weekData = retime(dayData,'weekly','sum');
weekEffort = retime(binEffort,'weekly','sum');
weekTable = retime(dayTable,'weekly','sum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create plots, binlog and pplog files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot Inter-Click Interval
iciIdx = find(ICIall > p.iciRange(1) & ICIall < p.iciRange(2));
% statistics
miciSel = mean(ICIall(iciIdx));
sdiciSel = std(ICIall(iciIdx));
meiciSel = median(ICIall(iciIdx));

figure(1); set(1,'name','Inter-Detection Interval')
h1 = gca;
centerIci = p.iciRange(1):1:p.iciRange(2);
[nhist] = histc(ICIall(iciIdx),centerIci);
bar(h1,centerIci,nhist, 'barwidth', 1, 'basevalue', 1)
xlim(h1,p.iciRange);
title(h1,sprintf('%s N=%d',p.filePrefix,length(ICIall)), 'Interpreter', 'none')
xlabel(h1,'Inter-Detection Interval (ms)')
ylabel(h1,'Counts')
% create labels and textbox
[~,modeIdx] = max(nhist);
moiciSel = centerIci(modeIdx);
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% save ici data and figure
icifn = [p.filePrefix,'_',p.speName,'_ici'];
saveas(h1,fullfile(p.tpwsDir,icifn),'fig')


%% Plot Peak-to-peak per click
% statistics
mpp = mean(PPall);
sdpp = std(PPall);
mepp = median(PPall);

% Plot histogram
figure(2); set(2,'name','Received Levels')
h2 = gca;
center = p.threshRL:1:p.rlHi;
[nhist] = histc(PPall,center);
bar(h2,center, nhist, 'barwidth', 1, 'basevalue', 1)
title(h2,sprintf('%s N=%d',p.filePrefix,length(PPall)), 'Interpreter', 'none')
xlabel(h2,'Peak-Peak Amplitude (dB)')
ylabel(h2,'Click Counts')
xlim(h2,[p.threshRL, p.rlHi])
% create labels and textbox
[~,modeIdx] = max(nhist);
mopp = center(modeIdx);
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
ppfn = [p.filePrefix,'_',p.speName,'_pp'];
save(fullfile(p.tpwsDir,ppfn),'clickData','PPall','center')
saveas(h2,fullfile(p.tpwsDir,ppfn),'fig')

% to store
centerHist = center;
histBins = nhist;
percHist = nhist*100/length(PPall);
%% Plot Peak-to-peak per bin
% statistics
mppBin = mean(binData.maxPP);
sdppBin = std(binData.maxPP);
meppBin = median(binData.maxPP);

% Plot histogram
figure(3); set(3,'name','Received Levels per Bin')
h3 = gca;
centerBin = p.threshRL:1:p.rlHi;
[nhistBin] = histc(binData.maxPP,centerBin);
bar(h3,centerBin, nhistBin, 'barwidth', 1, 'basevalue', 1)
title(h3,sprintf('%s N=%d',p.filePrefix,length(binData.maxPP)), 'Interpreter', 'none')
xlabel(h3,'Peak-Peak Amplitude (dB)')
ylabel(h3,[num2str(p.binDur),' min Bin Counts'])
xlim(h3,[p.threshRL, p.rlHi])
% create labels and textbox
[~,modeIdx] = max(nhist);
moppBin = centerBin(modeIdx);
mnlabelBin = sprintf('Mean = %0.2f', mppBin);
stdlabelBin = sprintf('Std = %0.2f', sdppBin);
melabelBin = sprintf('Median = %0.2f', meppBin);
molabelBin = sprintf('Mode = %0.2f', moppBin);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabelBin,stdlabelBin,...
    melabelBin,molabelBin});

% Save plot
binfn = [p.filePrefix,'_',p.speName,'_bin'];
save(fullfile(p.tpwsDir,binfn),'binData','centerBin')
saveas(h3,fullfile(p.tpwsDir,binfn),'fig')

% to store
percHistBin = nhistBin*100/length(binData.maxPP);
figure(5); set(5,'name','Weekly presence','DefaultAxesColor',[.8 .8 .8],'Position',[100 100 1400 600]) 
set(gca,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h5(1) = subplot(2,1,1);
h5(2) = subplot(2,1,2);
hold(h5(1), 'on')
bar(h5(1),weekTable.tbin,weekTable.Effort_Sec./604800*100,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(1),weekTable.tbin,weekTable.Count_Click./weekTable.Effort_Sec*100,'BarWidth', 1)
hold(h5(1), 'off')
hold(h5(2), 'on')
bar(h5(2),weekTable.tbin,weekTable.Effort_Bin./2016*100,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(2),weekTable.tbin,weekTable.Count_Bin./weekTable.Effort_Bin*100,'BarWidth', 1)
hold(h5(2), 'off')
ylabel(h5(1),{'Weekly Presence';'(% of clicks per week)'})
ylabel(h5(2),{'Weekly Presence';['(of 5 min bins per week)']})
title(h5(1),'Click Counting')
title(h5(2),'Group Counting');

ylim(h5(1),[0 round(max(weekTable.Count_Click./weekTable.Effort_Sec*100),2,'significant')*1.2])
ylim(h5(2),[0 round(max(weekTable.Count_Bin./weekTable.Effort_Bin*100),2,'significant')*1.2])
xlim(h5(1),[datetime(p.referenceTime),weekTable.tbin(end)+3])
xlim(h5(2),[datetime(p.referenceTime),weekTable.tbin(end)+3])
% define step according to number of weeks
if length(weekTable.tbin) > 53 && length(weekTable.tbin) <= 104 % 2 years
    step = calmonths(1);
elseif length(weekTable.tbin) > 104 && length(weekTable.tbin) <= 209 % 4 years
    step = calmonths(2);
elseif length(weekTable.tbin) > 209 && length(weekTable.tbin) <= 313 % 6 years
    step = calmonths(3);
elseif length(weekTable.tbin) > 313 && length(weekTable.tbin) <= 417 % 8 years
    step = calmonths(4);
elseif length(weekTable.tbin) >= 417 % more than 8 years
    step = calyears(1);
end
% define tick steps only if more than 1 year of data
if length(weekTable.tbin) > 53
    xtickformat(h5(1),'MMMyy')
    xtickformat(h5(2),'MMMyy')
    xticks(h5(1),[datetime(p.referenceTime),weekTable.tbin(1):step:weekTable.tbin(end)])
    xticks(h5(2),[datetime(p.referenceTime),weekTable.tbin(1):step:weekTable.tbin(end)])
    xtickangle(h5(1),45)
    xtickangle(h5(2),45)
end

weekfn = [p.filePrefix,'_',p.speName,'_weekly_presence'];
saveas(h5,fullfile(p.tpwsDir,weekfn),'fig')

summaryfn = [p.filePrefix,'_',p.speName,'_summaryData'];
save(fullfile(p.tpwsDir,summaryfn),'dayTable','weekTable','binData','NktTkt',...
    'NktTktbin','centerHist','percHist','percHistBin','histBins')

end

