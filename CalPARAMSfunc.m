function [ ] = CalPARAMSfunc(ct,cpp,csn,filePrefix,sp,sdir,detfn,p,tf,srate)
% Calculate inter-click interval, peak-to-peak, peak frequency, center 
% frequency, 3dB and 10 dB bandwidth, duration and rms amplitude.
% Create a figure of all the parameters.

%jah 5-6-14 modified 12-15-14 modified 5-23-15 for Kogia
% 8-7-15 modified for cpp 8-3-2016 modified for Simone tfParameters
% 1-16-2017 improved to work with detMOD
% 8-28-2017 modified by asb to work with grateful detEdit

%% Calculate teager energy and start end click
posClick = []; rmsClick = [];
for Cidx = 1:size(csn,1)
    hpdata = csn(Cidx,:).';
    teagerH = [];
    energy = spTeagerEnergy(hpdata);

    % Since we are operating on the high pass data, we'll
    % set the delay to zero.
    [SClicks, CClicks, noise, SNR] = dtHighResClick(srate*1000, energy, 0, hpdata, ...
                                                 teagerH, 0.001); % why 0.001? in triton is 0.015
    if isempty(CClicks)
        CClicks = [NaN NaN];
        SClicks = [NaN NaN];
    end
    posClick = [posClick;CClicks(1,:)];    
    rmsClick = [rmsClick;SClicks(1,:)]; 

    %plot teager energy and start and of single (SClicks) or full (CClicks)
    %click
    y = zeros(length(SClicks),1);

%     figure(100), plot(energy), hold on
%     plot(SClicks,y,'or')
%     if ~isnan(CClicks(1,1))
%         plot(CClicks,y,'og'), hold off
%     end
end
%Get rid of NaN data without duration
[nonan] = find(~isnan(posClick(:,1)) == 1);

%% Calculate parameters: rms, pf, cf, 3db, 10db, dur
[rms,peakFr,bw10db,bw3db,F0,dur] = ...
    paramDetEdit(csn,posClick,rmsClick,srate,tf,p.N);

%*****************************************************************
%% Peak-to-peak
% apply user range for db
if isempty(p.dbRange)
    p.dbRange = [min(cpp) max(cpp)];
end
% statistics
mpp = mean(cpp); 
sdpp = std(cpp);
mopp = mode(cpp); 
mepp = median(cpp);

% Plot histogram
figure(23);
nbinsdb = (p.dbRange(1):p.dbRange(2));
[y,centers] = hist(cpp,nbinsdb);
bar(centers,y)
title(sprintf('N=%d',length(cpp)))
xlabel('Peak-Peak Amplitude (dB)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
ppfn = strrep(detfn,'TPWS','pp');
saveas(fullfile(sdir,ppfn),'m') % save as matlab fig

%% RMS amplitude
% apply user range for db
if isempty(p.dbRange)
    p.dbRange = [min(rms) max(rms)];
end
% statistics
mrms = mean(rms); 
sdrms = std(rms);
morms = mode(rms); 
merms = median(rms);

% Plot histogram
figure(27);
nbinsdb = (p.dbRange(1):p.dbRange(2));
[y,centers] = hist(rms,nbinsdb);
bar(centers,y)
title(sprintf('N=%d',length(rms)))
xlabel('RMS Amplitude (dB)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mrms);
stdlabel = sprintf('Std = %0.2f', sdrms);
melabel = sprintf('Median = %0.2f', merms);
molabel = sprintf('Mode = %0.2f', morms);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
rmsfn = strrep(detfn,'TPWS','rms');
saveas(fullfile(sdir,rmsfn),'m')

%% Peak frequency
% statistics
mpeakFr = mean(peakFr(idur));
sdpeakFr = std(peakFr(idur));
mepeakFr = median(peakFr(idur));
mopeakFr = mode(peakFr(idur));

% Plot histogram
figure(24)
nbinsfr = (p.frRange(1):p.frRange(2));
[y,centers] = hist(peakFr(idur),nbinsfr);
bar(centers,y);
title(sprintf('N=%d',length(peakFr(idur))));
xlabel('Peak Frequency (kHz)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpeakFr);
stdlabel = sprintf('Std = %0.2f', sdpeakFr);
melabel = sprintf('Median = %0.2f', mepeakFr);
molabel = sprintf('Mode = %0.2f', mopeakFr);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
pffn = strrep(detfn,'TPWS','pf');
saveas(fullfile(sdir,pffn),'m')

%% Center frequency 
% statistics
mF0 = mean(F0(idur));
sdF0 = std(F0(idur));
meF0 = median(F0(idur));
moF0 = mode(F0(idur));

% Plot histogram
figure(28)
nbinsfr = (p.frRange(1):p.frRange(2));
[y,centers] = hist(F0(idur),nbinsfr);
bar(centers,y);
title(sprintf('N=%d',length(F0(idur))));
xlabel('Center Frequency (kHz)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mF0);
stdlabel = sprintf('Std = %0.2f', sdF0);
melabel = sprintf('Median = %0.2f', meF0);
molabel = sprintf('Mode = %0.2f', moF0);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
pffn = strrep(detfn,'TPWS','cf');
saveas(fullfile(sdir,pffn),'m')

%% 3 dB Bandwidth
% statistics
mbw3db = mean(bw3db(idur,3));
sdbw3db = std(bw3db(idur,3));
mebw3db = median(bw3db(idur,3));
mobw3db = mode(bw3db(idur,3));

% Plot histogram
figure(25)
nbinsdbw = (p.frdbwRange(1):p.frdbwRange(2));
[y,centers] = hist(bw3db(idur,3),nbinsdbw);
bar(centers,y);
title(sprintf('N=%d',length(bw3db(idur,3))));
xlabel('3 dB Bandwidth (kHz)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mbw3db);
stdlabel = sprintf('Std = %0.2f', sdbw3db);
melabel = sprintf('Median = %0.2f', mebw3db);
molabel = sprintf('Mode = %0.2f', mobw3db);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
bw3dbfn = strrep(detfn,'TPWS','3db');
saveas(fullfile(sdir,bw3dbfn),'m')

%% 10 dB Bandwidth
% statistics
mbw10db = mean(bw10db(idur,3));
sdbw10db = std(bw10db(idur,3));
mebw10db = median(bw10db(idur,3));
mobw10db = mode(bw10db(idur,3));

% Plot histogram
figure(29)
nbinsdbw = (p.frdbwRange(1):p.frdbwRange(2));
[y,centers] = hist(bw10db(idur,3),nbinsdbw);
bar(centers,y);
title(sprintf('N=%d',length(bw10db(idur,3))));
xlabel('10 dB Bandwidth (kHz)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mbw10db);
stdlabel = sprintf('Std = %0.2f', sdbw10db);
melabel = sprintf('Median = %0.2f', mebw10db);
molabel = sprintf('Mode = %0.2f', mobw10db);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
bw10dbfn = strrep(detfn,'TPWS','10db');
saveas(fullfile(sdir,bw10dbfn),'m')

%% Duration
% apply user range for db
if isempty(p.durRange)
    p.durRange = [min(dur) max(dur)];
end
% statistics
mdur = mean(dur(idur));
sddur = std(dur(idur));
medur = median(dur(idur));
modur = mode(dur(idur));

% Plot histogram
figure(26)
nbinsdur = (p.durRange(1):durstep:p.durRange(2));
[y,centers] = hist(dur(idur),nbinsdur);
bar(centers,y);
title(sprintf('N=%d',length(dur(idur))));
xlabel('Duration (\mus)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mdur);
stdlabel = sprintf('Std = %0.2f', sddur);
melabel = sprintf('Median = %0.2f', medur);
molabel = sprintf('Mode = %0.2f', modur);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
durfn = strrep(detfn,'TPWS','dur');
saveas(fullfile(sdir,durfn),'m')

%% Inter-click interval
ici = diff(ct)*24*60*60*1000; % in ms

% apply user range for ici
if isempty(p.iciRange)
    p.iciRange = [min(ici) max(ici)];
end
iciSel = ici(ici > p.iciRange(1) & ici < p.iciRange(2));

% statistics
miciSel = mean(iciSel);
sdiciSel = std(iciSel);
moiciSel = mode(iciSel);
meiciSel = median(iciSel);

% plot ici histogram
figure(22);
nbinsici = (p.iciRange(1):p.iciRange(2));
[y,centers] = hist(iciSel,nbinsici);
bar(centers,y);
xlim(p.iciRange);
title(sprintf('N=%d',length(ct)))
xlabel('Inter-Pulse Interval (ms)')
ylabel('Counts')
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% save ici data and figure
icifn = strrep(detfn,'TPWS','ici');
saveas(h22,fullfile(sdir,icifn),'m')

%% multiple plots
figure(30)
% peak-to-peak
subplot(2,4,1),hist(cpp,nbinsdb,'k')
xlabel('peak-peak amplitude dB')
xlim([p.dbRange(1),p.dbRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(cpp)))

% RMS amplitude
subplot(2,4,5),hist(rms,p.dbRange)
xlabel('RMS amplitude dB')
xlim([p.dbRange(1),p.dbRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(rms)))

% peak frequency
subplot(2,4,2),hist(peakFr(idur),nbinsfr)
xlabel('peak frequency (kHz)')
xlim([p.frRange(1),p.frRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(peakFr(idur))));

% center frequency
subplot(2,4,6),hist(F0(idur),nbinsfr)
xlabel('peak frequency (kHz)')
xlim([p.frRange(1),p.frRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(F0(idur))));

% 3db bandwidth
subplot(2,4,3),hist(bw3db(idur,3),nbinsdbw)
xlabel('-3dB Bandwidth(kHz)')
xlim([p.frdbwRange(1),p.frdbwRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(bw3db(idur,3))));

% 10db bandwidth
subplot(2,4,7),hist(bw10db(idur,3),nbinsdbw)
xlabel('-10dB Bandwidth(kHz)')
xlim([p.frdbwRange(1),p.frdbwRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(bw10db(idur,3))));

% duration
subplot(2,4,4),hist(dur(idur),nbinsdur)
xlabel('duration (\mus)')
xlim([p.durRange(1),p.durRange(2)]);
ylabel('counts')
title(sprintf('N=%d',length(dur(idur))));

% inter-click interval
subplot(2,4,8),hist(iciSel,nbinsici)
xlim(p.iciRange);
xlabel('inter-click interval (ms)')
ylabel('counts')
title(sprintf('N=%d',length(ct)))

% save figure
allfn = strrep(detfn,'TPWS','params');
save(fullfile(sdir,allfn),'rms','peakFr','F0','bw3db','bw10db','dur','iciSel')
saveas(fullfile(sdir,allfn),'m')
saveas(fullfile(sdir,allfn),'png') % also save in png

%*****************************************************************
% Create excel with data parameters
vecTimes = datevec(ct(1:end));
ici(end+1) = ici(end);
xlsParams = [vecTimes, 1000*ici',cpp,rms,peakFr,F0,bw3db,bw10db,dur'];
%xlswrite(fn3,xlsParams);  % write time and click count by bin data to XLS
% Open Excel, add workbook, change active worksheet,
% get/put array, save, and close
% First open an Excel Server
Excel = actxserver('Excel.Application');
set(Excel, 'Visible', 1);
% Insert a new workbook
Workbooks = Excel.Workbooks;
Workbook = invoke(Workbooks, 'Add');
% Make the second sheet active
Sheets = Excel.ActiveWorkBook.Sheets;
sheet2 = get(Sheets, 'Item', 1);
invoke(sheet2, 'Activate');
% Get a handle to the active sheet
Activesheet = Excel.Activesheet;
% Put a MATLAB array into Excel
%A = [1 2; 3 4];  
A = xlsParams;
% B = ['Year','Month','Day','Hour','Min','Sec','ici(ms)',...
%     'cpp(dB)','peakFr(kHz)','centerFr(kHz)','bw3db(kHz)','bw10db(kHz)',...
%     'dur(ns)'];
la = length(A);
rstrg = ['A1:N',num2str(la)];
ActivesheetRange = get(Activesheet,'Range',rstrg);
set(ActivesheetRange, 'Value', A);
% Get back a range.  It will be a cell array, 
% since the cell range can
% contain different types of data.
% Range = get(Activesheet, 'Range', 'A1:B2');
% B = Range.value;
% % Convert to a double matrix.  The cell array must contain only scalars.
% B = reshape([B{:}], size(B));
% Save the workbook
xlsfn = strrep(allfn,'.mat','.xls');
invoke(Workbook, 'SaveAs', xlsfn);
% To avoid saving the workbook and being prompted to do so,
% uncomment the following code.
Workbook.Saved = 1;
invoke(Workbook, 'Close');
% Quit Excel
invoke(Excel, 'Quit');
% End process
delete(Excel);

