function SumPPICIBin(varargin)

%% get user input and set up file names
n = 1;
while n <= length(varargin)
    switch varargin{n}
        case 'filePrefix'
            filePrefix = varargin{n+1}; n=n+2;
        case 'sp'
            sp = varargin{n+1}; n=n+2;
        case 'sdir'
            sdir = varargin{n+1}; n=n+2;
        case 'concatFiles'
            concatFiles = varargin{n+1}; n=n+2;
        case 'effort'
            effort = varargin{n+1}; n=n+2;
        case 'referenceTime'
            refTime = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end
sdir = [sdir,'\SumPPICI'];
mkdir(sdir)  % create place to save files
% cd([tpwsPath,'\SumPPICI'])
%% get default parameters
p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin');

%% Concatenate variables
PPall = []; TTall = []; ICIall = []; % initialize matrices
parfor idsk = 1 : length(concatFiles)
    % Load file
    fprintf('Loading %d/%d file %s\n',idsk,length(concatFiles),char(concatFiles{idsk}))
    D = load(char(concatFiles{idsk}));
    
    % find times outside effort (sometimes there are detections
    % which are from the audio test at the beggining of the wav file)
    within = cell2mat(arrayfun(@(x)sum(isbetween(x,datenum(effort.Start),datenum(effort.End))),D.MTT,'uni',false));
    goodIdx = find(within ~= 0);
    MTT = D.MTT(goodIdx); % only keep the good detections
    MPP = D.MPP(goodIdx);
    if isrow(MTT)
        MTT = MTT';
    end
    if isrow(MPP)
        MPP = MPP';
    end
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

secMonitEffort = sum(effort.diffSec);
binMonitEffort = sum(effort.roundbin);

% convert intervals in bins 
binEffort = intervalToBinTimetable(effort.Start,effort.End,p); 
binEffort.sec = binEffort.bin*(p.binDur*60);

%% get average of detection by effort
NktTkt = positiveCounts/secMonitEffort;
NktTktbin = positiveBins/binMonitEffort;
disp(['Nkt/Tkt for click counting: ',num2str(NktTkt)]);
disp(['Nkt/Tkt for group counting: ',num2str(NktTktbin)]);

%% group data by days and weeks
Click = retime(binData(:,1),'daily','sum'); % #click per day
Bin = retime(binData(:,1),'daily','count'); % #bin per day

dayData = synchronize(Click,Bin);
dayEffort = retime(binEffort,'daily','sum');
dayTable = synchronize(dayData,dayEffort);
dayTable.Properties.VariableNames{'bin'} = 'Effort_Bin';
dayTable.Properties.VariableNames{'sec'} = 'Effort_Sec';

weekData = retime(dayData,'weekly','mean');
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
moiciSel = mode(ICIall(iciIdx));
meiciSel = median(ICIall(iciIdx));

figure(1); set(1,'name','Inter-Click Interval')
h1 = gca;
centerIci = p.iciRange(1):1:p.iciRange(2);
[nhist] = histc(ICIall(iciIdx),centerIci);
bar(h1,centerIci,nhist, 'barwidth', 1, 'basevalue', 1)
xlim(h1,p.iciRange);
title(h1,sprintf('%s N=%d',filePrefix,length(ICIall)), 'Interpreter', 'none')
xlabel(h1,'Inter-Click Interval (ms)')
ylabel(h1,'Counts')
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% save ici data and figure
icifn = [filePrefix,'_',p.speName,'_ici'];
saveas(h1,fullfile(sdir,icifn),'fig')


%% Plot Peak-to-peak per click
% statistics
mpp = mean(PPall);
sdpp = std(PPall);
mopp = mode(PPall);
mepp = median(PPall);

% Plot histogram
figure(2); set(2,'name','Received Levels')
h2 = gca;
center = p.threshRL:1:p.p1Hi;
[nhist] = histc(PPall,center);
bar(h2,center, nhist, 'barwidth', 1, 'basevalue', 1)
title(h2,sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
xlabel(h2,'Peak-Peak Amplitude (dB)')
ylabel(h2,'Click Counts')
xlim(h2,[p.threshRL, p.p1Hi])
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
ppfn = [filePrefix,'_',p.speName,'_pp'];
save(fullfile(sdir,ppfn),'clickData','PPall','center')
saveas(h2,fullfile(sdir,ppfn),'fig')

%% Plot Peak-to-peak per bin
% statistics
mppBin = mean(binData.maxPP);
sdppBin = std(binData.maxPP);
moppBin = mode(binData.maxPP);
meppBin = median(binData.maxPP);

% Plot histogram
figure(3); set(3,'name','Received Levels per Bin')
h3 = gca;
centerBin = p.threshRL:1:p.p1Hi;
[nhistBin] = histc(binData.maxPP,centerBin);
bar(h3,centerBin, nhistBin, 'barwidth', 1, 'basevalue', 1)
title(h3,sprintf('%s N=%d',filePrefix,length(binData.maxPP)), 'Interpreter', 'none')
xlabel(h3,'Peak-Peak Amplitude (dB)')
ylabel(h3,[num2str(p.binDur),' min Bin Counts'])
xlim(h3,[p.threshRL, p.p1Hi])
% create labels and textbox
mnlabelBin = sprintf('Mean = %0.2f', mppBin);
stdlabelBin = sprintf('Std = %0.2f', sdppBin);
melabelBin = sprintf('Median = %0.2f', meppBin);
molabelBin = sprintf('Mode = %0.2f', moppBin);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabelBin,stdlabelBin,...
    melabelBin,molabelBin});

% Save plot
binfn = [filePrefix,'_',p.speName,'_bin'];
save(fullfile(sdir,binfn),'binData','centerBin')
saveas(h3,fullfile(sdir,binfn),'fig')

%% Plot weekly mean of detections 
figure(5); set(5,'name','Weekly presence','DefaultAxesColor',[.8 .8 .8],'Position',[100 100 1400 600]) 
set(gca,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h5(1) = subplot(2,1,1);
h5(2) = subplot(2,1,2);
hold(h5(1), 'on')
bar(h5(1),weekEffort.tbin,weekEffort.sec,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(1),weekData.tbin,weekData.Count_Click,'BarWidth', 1)
hold(h5(1), 'off')
hold(h5(2), 'on')
bar(h5(2),weekEffort.tbin,weekEffort.bin,'FaceColor',[1 1 1],'EdgeColor','none','BarWidth', 1)
bar(h5(2),weekData.tbin,weekData.Count_Bin,'BarWidth', 1)
hold(h5(2), 'off')
%set(h5(1),'xticklabel', '');
ylabel(h5(1),{'Weekly mean';'of clicks per day'})
ylabel(h5(2),{'Weekly mean';['of ', num2str(p.binDur), ' min bins per day']})
title(h5(1),'Click Counting')
title(h5(2),'Group Counting');
% axis (h5(1),'tight')
% axis (h5(2),'tight')

ylim(h5(1),[0 round(max(weekData.Count_Click),2,'significant')*1.2])
ylim(h5(2),[0 round(max(weekData.Count_Bin),2,'significant')*1.2])
xlim(h5(1),[datetime(refTime),weekEffort.tbin(end)+31])
xlim(h5(2),[datetime(refTime),weekEffort.tbin(end)+31])
% define step according to number of weeks
if length(weekData.tbin) > 53 && length(weekData.tbin) <= 104 % 2 years
    step = calmonths(1);
elseif length(weekData.tbin) > 104 && length(weekData.tbin) <= 209 % 4 years
    step = calmonths(2);
elseif length(weekData.tbin) > 209 && length(weekData.tbin) <= 313 % 6 years
    step = calmonths(3);
elseif length(weekData.tbin) > 313 && length(weekData.tbin) <= 417 % 8 years
    step = calmonths(4);
elseif length(weekData.tbin) >= 417 % more than 8 years
    step = calyears(1);
end
% define tick steps only if more than 1 year of data
if length(weekData.tbin) > 53
    xtickformat(h5(1),'MMMyy')
    xtickformat(h5(2),'MMMyy')
    xticks(h5(1),datetime(refTime):step:weekEffort.tbin(end)+31)
    xticks(h5(2),datetime(refTime):step:weekEffort.tbin(end)+31)
    xtickangle(h5(1),45)
    xtickangle(h5(2),45)
end

weekfn = [filePrefix,'_',p.speName,'_weekly_presence'];
saveas(h5,fullfile(sdir,weekfn),'fig')

summaryfn = [filePrefix,'_',p.speName,'_summaryData_forDensity'];
save(fullfile(sdir,summaryfn),'dayTable','weekTable','refTime','NktTkt','NktTktbin')

end

