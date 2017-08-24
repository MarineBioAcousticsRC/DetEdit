function [ ] = CalPARAMSfunc(itnumo,MTT,MPP,MSN,stn,dpn,spe,...
    pn1,fs,tffreq,tfuppc)
% Calculate the ICI and the PP 
%jah 5-6-14 modified 12-15-14 modified 5-23-15 for Kogia
% 8-7-15 modified for pp 8-3-2016 modified for Simone tfParameters
% 1-16-2017 improved to work with detMOD
% DEFINE FILES
% ici=inter-click-interval; pp = peak-peak
icifn = [stn,dpn,'_',spe,'_ici',num2str(itnumo)];
icifnfig = [stn,dpn,'_',spe,'_ici',num2str(itnumo),'.pdf'];
xlsfn = [stn,dpn,'_',spe,'params.xls'];
ppfn = [stn,dpn,'_',spe,'_pp',num2str(itnumo)];
ppfnfig = [stn,dpn,'_',spe,'_pp',num2str(itnumo),'.pdf'];
% pf = peak-frequency; cf = center-frequency
pffn = [stn,dpn,'_',spe,'_pf',num2str(itnumo)];
pffnfig = [stn,dpn,'_',spe,'_pf',num2str(itnumo),'.pdf'];
cffn = [stn,dpn,'_',spe,'_cf',num2str(itnumo)];
cffnfig = [stn,dpn,'_',spe,'_cf',num2str(itnumo),'.pdf'];
% 3db = 3dB-bandwidth; 10db = 10dB-bandwidth
bw3dbfn = [stn,dpn,'_',spe,'_3db',num2str(itnumo)];
bw3dbfnfig = [stn,dpn,'_',spe,'_3db',num2str(itnumo),'.pdf'];
bw10dbfn = [stn,dpn,'_',spe,'_10db',num2str(itnumo)];
bw10dbfnfig = [stn,dpn,'_',spe,'_10db',num2str(itnumo),'.pdf'];
% dur = duration  rms= rms ampltidue
durfn = [stn,dpn,'_',spe,'_dur',num2str(itnumo)];
durfnfig = [stn,dpn,'_',spe,'_dur',num2str(itnumo),'.pdf'];
rmsfn = [stn,dpn,'_',spe,'_rms',num2str(itnumo)];
rmsfnfig = [stn,dpn,'_',spe,'_rms',num2str(itnumo),'.pdf'];
% all plots
allfn = [stn,dpn,'_',spe,'_all',num2str(itnumo)];
allfnfig = [stn,dpn,'_',spe,'_all',num2str(itnumo),'.pdf'];
% fn = fullfile(detpn,detfn);
fn1 = fullfile(pn1,'\',icifn);
fn2 = fullfile(pn1,'\',icifnfig);
fn3 = fullfile(pn1,'\',xlsfn);
fn4 = fullfile(pn1,'\',ppfn);
fn5 = fullfile(pn1,'\',ppfnfig);
fn6 = fullfile(pn1,'\',pffn);
fn7 = fullfile(pn1,'\',pffnfig);
fn8 = fullfile(pn1,'\',cffn);
fn9 = fullfile(pn1,'\',cffnfig);
fn10 = fullfile(pn1,'\',bw3dbfn);
fn11 = fullfile(pn1,'\',bw3dbfnfig);
fn12 = fullfile(pn1,'\',bw10dbfn);
fn13 = fullfile(pn1,'\',bw10dbfnfig);
fn14 = fullfile(pn1,'\',durfn);
fn15 = fullfile(pn1,'\',durfnfig);
fn16 = fullfile(pn1,'\',rmsfn);
fn17 = fullfile(pn1,'\',rmsfnfig);
ct = MTT; 
pp = MPP; % needs to be TF corrected version of pp
sn = MSN;
%Calculate duration
posClick = []; rmsClick = [];
for Cidx = 1:size(MSN,1)
    hpdata = MSN(Cidx,:).';
    teagerH = [];
    energy = spTeagerEnergy(hpdata);

    % Since we are operating on the high pass data, we'll
    % set the delay to zero.
    [SClicks, CClicks, noise, SNR] = dtHighResClick(fs, energy, 0, hpdata, ...
                                                 teagerH, 0.001);
    if isempty(CClicks)
        CClicks = [NaN NaN];
        SClicks = [NaN NaN];
    end
    posClick = [posClick;CClicks(1,:)];    
%      display([num2str(CClicks(1)),' click #',num2str(Cidx),...
%          '  ',num2str(size(posClick))]);
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
%
% CALCULTE rms pf cf 3db 10db dur
[rms,peakFr,bw10db,bw3db,F0,dur] = ...
    paramDetEdit(MSN,posClick,rmsClick,fs,tffreq,tfuppc); 
idur = find(~isnan(dur) > 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot histograms in one Figure(30)
figure(30)
if strcmp(spe, 'Kogia')
    fpcrange = 80: 1 : 160 ; % freq range for peak and center in kHz
    f310range = 0: 1 : 80 ; % freq range for -3 and -10 bw in kHz
    dbrange = 100 : 1 : 150; % dB range for pp and rms
    drange = 30 : 3 : 111; % duration range in us
    ipirange = 40 : 1 :130; % ici range in ms
end
if (strcmp(spe, 'Gervais') || strcmp(spe, 'Cuviers'))
    fpcrange = 0: 1 : 100 ; % freq range for peak and center in kHz
    f310range = 0: 1 : 80 ; % freq range for -3 and -10 bw in kHz
    dbrange = 90 : 1 : 170; % dB range for pp and rms
    drange = 30 : 2 : 300; % duration range in us
    ipirange = 40 : 1 :750; % ici range in ms
end
if (strcmp(spe, 'Beluga'))
    fpcrange = 0: 1 : 100 ; % freq range for peak and center in kHz
    f310range = 0: 1 : 80 ; % freq range for -3 and -10 bw in kHz
    dbrange = 90 : 1 : 170; % dB range for pp and rms
    drange = 10 : 2 : 300; % duration range in us
    ipirange = 20 : 1 :500; % ici range in ms
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Save peak-peak data
figure(23)
[y,centers] = hist(pp,dbrange);
bar(centers,y)
mpp = mean(pp); 
sdpp = std(pp);
imax = find(y == max(y)); mepp = median(pp);
ppstr=['Mean= ',num2str(mpp),' Median= ',num2str(mepp),...
    ' Mode= ',num2str(dbrange(imax)),...
    ' StDev= ',num2str(sdpp),' Number= ',num2str(length(pp))];
title(ppstr)
xlabel('Peak-Peak Amplitude (dB)')
save(fn4,'pp')
saveas(gcf,fn5,'pdf')
figure(30)
subplot(2,4,1),hist(pp,dbrange,'k')
xlabel('peak-peak amplitude dB')
xlim([min(dbrange),max(dbrange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(pp)));
title(strp);

% Save pf data
figure(24)
[y,centers] = hist(peakFr(idur),fpcrange);
bar(centers,y);
mpeakFr = mean(peakFr(idur));
sdpeakFr = std(peakFr(idur));
imax = find(y == max(y)); mepeakFr = median(peakFr(idur));
ppstr=['Mean= ',num2str(mpeakFr),'Median= ',num2str(mepeakFr),...
    'Mode= ',fpcrange(imax),...
    'dB  StDev= ',num2str(sdpeakFr),' Number= ',...
    num2str(length(peakFr(idur)))];
title(ppstr)
xlabel('Peak Frequency (kHz)')
save(fn6,'peakFr')
saveas(gcf,fn7,'pdf')
figure(30)
subplot(2,4,2),hist(peakFr(idur),fpcrange)
xlabel('peak frequency (kHz)')
xlim([min(fpcrange),max(fpcrange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(peakFr(idur))));
title(strp);

% Save 3 dB bandwidth data
figure(25)
[y,centers] = hist(bw3db(idur,3),f310range);
bar(centers,y);
mbw3db = mean(bw3db(idur,3));
sdbw3db = std(bw3db(idur,3));
imax = find(y == max(y)); mebw3 = median(bw3db(idur,3));
ppstr=[' Mean= ',num2str(mbw3db),' Median= ',num2str(mebw3),...
    ' Mode= ',num2str(f310range(imax)),...
    'dB  StDev= ',num2str(sdbw3db),' Number= ',...
    num2str(length(bw3db(idur,3)))];
title(ppstr)
xlabel('3 dB Bandwidth (kHz)')
save(fn10,'bw3db')
saveas(gcf,fn11,'pdf')
figure(30)
subplot(2,4,3),hist(bw3db(idur,3),f310range)
xlabel('-3dB Bandwidth(kHz)')
xlim([min(f310range),max(f310range)]);
ylabel('counts')
strp=strcat('n=',num2str(length(bw3db(idur,3))));
title(strp);

% Save dur data
figure(26)
[y,centers] = hist(dur(idur),drange);
bar(centers,y);
mdur = mean(dur(idur));
sddur = std(dur(idur));
imax = find(y == max(y)); medur = median(dur(idur));
ppstr=[' Mean= ',num2str(mdur),' Median= ',num2str(medur),...
    ' Mode= ',num2str(drange(imax)),...
    'dB  StDev= ',num2str(sddur),' Number= ',...
    num2str(length(dur(idur)))];
title(ppstr)
xlabel('Duration(ns)')
save(fn14,'dur')
saveas(gcf,fn15,'pdf')
figure(30)
subplot(2,4,4),hist(dur(idur),drange)
xlabel('duration (us)')
xlim([min(drange),max(drange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(dur(idur))));
title(strp);

% RMS amplitude
figure(27)
[y,centers] = hist(rms(idur),dbrange);
bar(centers,y);
mrms = mean(rms(idur));
sdrms = std(rms(idur));
imax = find(y == max(y)); merms = median(rms(idur));
ppstr=[' Mean= ',num2str(mrms),' Median= ',num2str(merms),...
    ' Mode= ',num2str(dbrange(imax)),...
    'dB  StDev= ',num2str(sdrms),' Number= ',num2str(length(rms(idur)))];
title(ppstr)
xlabel('RMS Amplitude (dB)')
save(fn16,'rms')
saveas(gcf,fn17,'pdf')
figure(30)
subplot(2,4,5),hist(rms,dbrange)
xlabel('RMS amplitude dB')
xlim([min(dbrange),max(dbrange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(rms)));
title(strp);

% Save cf data
figure(28)
[y,centers] = hist(F0(idur),fpcrange);
bar(centers,y);
mF0 = mean(F0(idur));
sdF0 = std(F0(idur));
imax = find(y == max(y)); meF0 = median(F0(idur));
ppstr=[' Mean= ',num2str(mF0),' Median= ',num2str(meF0),...
    ' Mode= ',num2str(fpcrange(imax)),...
    'dB  StDev= ',num2str(sdF0),' Number= ',...
    num2str(length(F0(idur)))];
title(ppstr)
xlabel('Center Frequency (kHz)')
save(fn8,'F0')
saveas(gcf,fn9,'pdf')
figure(30)
subplot(2,4,6),hist(F0(idur),fpcrange)
xlabel('Center Frequency (kHz)')
xlim([min(fpcrange),max(fpcrange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(F0(idur))));
title(strp);

% Save 10 dB bandwidth data
figure(29)
[y,centers] = hist(bw10db(idur,3),f310range);
bar(centers,y);
mbw10db = mean(bw10db(idur,3));
sdbw10db = std(bw10db(idur,3));
imax = find(y == max(y)); mebw10 = median(bw10db(idur,3));
ppstr=[' Mean= ',num2str(mbw10db),' Median= ',num2str(mebw10),...
    ' Mode= ',num2str(f310range(imax)),...
    'dB  StDev= ',num2str(sdbw10db),' Number= ',...
    num2str(length(bw10db(idur,3)))];
title(ppstr)
xlabel('10 dB Bandwidth (kHz)')
save(fn12,'bw10db')
saveas(gcf,fn13,'pdf')
figure(30)
subplot(2,4,7),hist(bw10db(idur,3),f310range)
xlabel('-10dB bandwidth (kHz)')
xlim([min(f310range),max(f310range)]);
ylabel('counts')
strp=strcat('n=',num2str(length(bw10db(idur,3))));
title(strp);

% CALCULATE and SAVE ICI
ici = diff(ct)*24*60*60*1000; % in ms
%define threshold for short ici to be 10 ms and long to be 300 ms
iciIdx = find(ici > 10 & ici < max(ipirange));
figure(22)
icif = ici(iciIdx);  % select in range
[y,centers] = hist(icif,ipirange);
bar(centers,y);
xlim([min(ipirange)+10 max(ipirange)-50]);
micif = mean(icif);
sdicif = std(icif);
imax = find(y == max(y)); meicif = median(icif);
ppstr=[' Mean= ',num2str(micif),' Median= ',num2str(meicif),...
    ' Mode= ',num2str(ipirange(imax)),...
    ' StDev= ',num2str(sdicif),' Number= ',num2str(length(icif))];
title(ppstr)
xlabel('Inter-Pulse Interval (ms)')
save(fn1,'icif')
saveas(gcf,fn2,'pdf')
figure(30)
subplot(2,4,8),hist(icif,ipirange)
xlim([min(ipirange)+10 max(ipirange)-50]);
% axis([min(ipirange)+10 max(ipirange) 0 1.1*max(y)]);
xlabel('inter-click interval (ms)')
xlim([min(ipirange),max(ipirange)]);
ylabel('counts')
strp=strcat('n=',num2str(length(icif)));
title(strp);
%*****************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xls for Danielle
icixls = datevec(ct(1:end));
ici(end+1) = ici(end);
icixls = [icixls, 1000*ici',pp,rms,peakFr,F0,bw3db,bw10db,dur'];
%xlswrite(fn3,icixls);  % write time and click count by bin data to XLS
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
A = icixls;
% B = ['Year','Month','Day','Hour','Min','Sec','ici(ms)',...
%     'pp(dB)','peakFr(kHz)','centerFr(kHz)','bw3db(kHz)','bw10db(kHz)',...
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
% Now save the workbook
%invoke(Workbook, 'SaveAs', 'myfile.xls');
invoke(Workbook, 'SaveAs', fn3);
% To avoid saving the workbook and being prompted to do so,
% uncomment the following code.
Workbook.Saved = 1;
invoke(Workbook, 'Close');
% Quit Excel
invoke(Excel, 'Quit');
% End process
delete(Excel);

