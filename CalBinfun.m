function [] = CalBinfun(tjan,t,dpnspe,pn1,binDur,gt,ise,seff,eeff)
% Changed a function with variable bin length 1-15-16
%
%dn = sscanf(dpn,['%2i']);  % make sure to get only deploy number
if (str2double(ise) > 0)
    dayfn = [dpnspe,ise,'_day.txt'];
    efffn = [dpnspe,ise,'_eff.txt'];
    binfn = [dpnspe,ise,'_bin.xls'];
else
    dayfn = [dpnspe,'_day.txt'];
    efffn = [dpnspe,'_eff.txt'];
    binfn = [dpnspe,'_bin.xls'];
end
fn1 = fullfile(pn1,'\',dayfn);
fn2 = fullfile(pn1,'\',efffn);
fn3 = fullfile(pn1,'\',binfn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find edges (start and end times) of bouts or sessions
dt = diff(t)*24*60*60; % time between detections 
%                           convert from days to seconds
I = [];
I = find(dt>gt);  % find start of gaps , gt in sec
sb = [t(1),t(I+1)];   % start time of bout
eb = [t(I),t(end)];   % end time of bout
dd = t(end)-t(1);     % deployment duration [d]
nb = length(sb);        % number of bouts
bd = (eb - sb)*24*60;      % duration of bout in min
figure(24)
hist(bd,100)
ylabel('Number of Bouts');
xlabel('Duration of Bouts in Minutes');
%
nd = length(t); %  nd = number of detections
% data for start and end of effort
% convert to MATLAB datenum
ts = x2mdate(seff);
te = x2mdate(eeff);
dur = te - ts;
% count number of time bins in data
% tjan = datenum([2009 1 1 0 0 0]);
tfirst = floor(t(1) - tjan);    %day of first detection
tlast = floor(t(end) - tjan);   %day of last detection
nday = tlast - tfirst + 1;          %num days with detections
tday = [tfirst:tlast]';             %tday array of #bins/day
tday(:,2) = 0; tday(:,3) = 0; % 2nd col is #bin, 3rd col is # clicks
% sort detections into time bins
%binDur      % bin duration [minutes]
nbin = ceil(dur*24*60/binDur); % num bins in entire deployment
nbind = ceil(24*60/binDur);     %num bins in one day
bin = 1;
kb = 1;
C = zeros(nbin,1); % the number of detections in time bin
T = zeros(nbin,1); % time for time bin
while kb <= nd % nd = number detections
    tv = datevec(t(kb));     % put time (kb) in vector format
    tbin = floor(tv(5)/binDur);  % define min time bin
    t0 = datenum([tv(1:4) tbin*binDur 00]); % lower bound
    t1 = datenum([tv(1:4) (tbin+1)*binDur 00]); % upper bound
    I = [];     % index for times in time bin
    I = find(t >= t0 & t < t1);
    if ~isempty(I)
        C(bin) = length(I); % the number of detections in time bin
        T(bin) = t0;        % time for time bin beginning
        iday = floor(t0 - tjan + 1 - tfirst);
        tday(iday,2) = tday(iday,2) + 1;  % count bins per day
        tday(iday,3) = tday(iday,3) + length(I);  % count clicks per day
        kb = I(end) + 1;      % set loop index to next time
        bin = bin + 1;  % increment bin number
    else
        if t(kb) == 0
            kb = kb + 1;
            bin = bin+1;
        end
    end
end
% filter empty and low number bins
minNdet = 1;        % minimum number of detections per bin
% if CHANGE to NOT = 1 ADD code below
KB = [];
KB = find(C >= minNdet);
if isempty(KB)
    disp(['No bins with at least ',num2str(minNdet),' detections'])
    binT = [0];
    binC = [0];
else
    %binT = T(KB) + datenum([0 0 0 0 binDur/2 0]); % move to middle of bin
    binT = T(KB);  % leave at beginning of bin
    binC = C(KB);
end
binxls = datevec(binT);
binxls = [binxls, binC];
% write time and click count by bin data to XLS
%xlswrite(fn3,binxls);
% Open Excel, add workbook, change active worksheet,
Excel = actxserver('Excel.Application');
set(Excel, 'Visible', 1);
Workbooks = Excel.Workbooks;
Workbook = invoke(Workbooks, 'Add');
Sheets = Excel.ActiveWorkBook.Sheets;
sheet2 = get(Sheets, 'Item', 1);
invoke(sheet2, 'Activate');
Activesheet = Excel.Activesheet;
A = binxls;
la = length(A);
rstrg = ['A1:G',num2str(la)];
ActivesheetRange = get(Activesheet,'Range',rstrg);
set(ActivesheetRange, 'Value', A);
invoke(Workbook, 'SaveAs', fn3);
% To avoid saving the workbook and being prompted to do so,
% uncomment the following code.
Workbook.Saved = 1;
invoke(Workbook, 'Close');
% Quit Excel
invoke(Excel, 'Quit');
% End process
delete(Excel);
% ADD THIS CODE BACK WHEN minNdet is NOT = 1
% ALSO SUBTRACT Number click tday(Bday,3)
% filter low click number bins from day counts
BB = find(C < minNdet & C > 0);
if ~isempty(BB)
    for iBB = 1:length(BB)
        Bday = floor(T(BB(iBB)) - tjan + 1 - tfirst);
        tday(Bday,2) = tday(Bday,2)  - 1;
        tday(Bday,3) = tday(Bday,3) - C(BB(iBB));
    end
    nClick = sum(C(BB));
    disp([' Removed: ',num2str(nClick),' Clicks with less than ',...
        num2str(minNdet),' per bin']);
end
% Compensate for partial days at beginning or end
% add zeros at beginning and end
sday = (ts - tjan);    %day of first effort
eday = (te - tjan);   %day of last effort
fsday = floor(sday);
feday = floor(eday);
sfrac = 1 - (sday - fsday);
efrac = (eday - feday);
if (tfirst == fsday)
    %removed since want #bins not #bin/day
%     disp([' Correct first day ',num2str(sfrac)]);
%     if (sfrac > 0)
%         tday(1,2) = round(tday(1,2)/sfrac);
%     end
else
    nzero = zeros(tfirst - fsday , 2);
    dadd = [fsday : (tfirst -1)]';
    tday = [dadd, nzero  ; tday];
end
nzero = [];
dadd = [];
if (tlast == feday)
     %removed since want #bins not #bin/day
%     disp([' Correct last day ',num2str(efrac)]);
%     if (sfrac > 0)
%         tday(length(tday),2) = round(tday(length(tday),2)/efrac);
%     end
else
    nzero = zeros(feday - tlast , 2);
    dadd = [(tlast +1) : feday]';
    tday = [tday; dadd, nzero];
end
% add effort as 4th column
efft = nbind*ones(feday - (fsday-1) , 1);
efft(1) = sfrac * efft(1);
efft(end) = efrac * efft(end);
if (ise >= 1) % needed when broken effort
    ifs = find(tday(:,1) == fsday);
    ife = find(tday(:,1) == feday);
    tdayefft = tday(ifs:ife,:);
    tday = tdayefft;
end
tday = [tday,efft];
%
nztday = find(tday(:,2) > 0);
disp(['Number of clicks = ',num2str(sum(binC))]);
disp(['Number of ',num2str(binDur),' min Bins Effort = ',...
    num2str(nbin)]);
save(fn2,'nbin','-ascii')
disp(['Number of Bins w/ ',num2str(minNdet),' clicks = ',...
    num2str(length(binT))]);
disp(['Number of Days with clicks  ',   num2str(length(nztday))]);
save(fn1,'tday','-ascii')
figure(10)
bar(tday(:,1),tday(:,3)/tday(:,2),'k') %average daily clicks/bin
title('Average of Daily clicks per bin')
figure(11)
plot(tday(:,3),tday(:,2),'ok') %num clicks vs num bin
ylabel('Number bin per day');
xlabel('Number Clicks per day');
%
disp(['This is the end']);
end

