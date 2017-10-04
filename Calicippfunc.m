function [ ] = Calicippfunc(MTT,MPP,MSP,sdir,detfn,p,srate)
% Calculate the ICI and the PP 
close all

%% Inter-Click Interval
ici = diff(MTT)*24*60*60*1000; % in ms

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
h22= figure(22);
nbinsici = (p.iciRange(1):p.iciRange(2));
[y,centers] = hist(iciSel,nbinsici);
bar(centers,y);
xlim(p.iciRange);
title(sprintf('N=%d',length(MTT)))
xlabel('Inter-Pulse Interval (ms)')
ylabel('Counts')
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});
axis tight

% save ici data and figure
icifn = strrep(detfn(1:end-4),'TPWS','ici');
%saveas(h22,fullfile(sdir,icifn))
saveas(h22,fullfile(sdir,icifn),'png')

%% Peak-to-peak
% apply user range for db
if isempty(p.dbRange)
    p.dbRange = [min(MPP) max(MPP)];
end
% statistics
mpp = mean(MPP); 
sdpp = std(MPP);
mopp = mode(MPP); 
mepp = median(MPP);

% Plot histogram
h23 = figure(23);
nbinsdb = (p.dbRange(1):p.dbRange(2));
[y,centers] = hist(MPP,nbinsdb);
bar(centers,y)
title(sprintf('N=%d',length(MPP)))
xlabel('Peak-Peak Amplitude (dB)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});
axis tight

% Save plot
ppfn = strrep(detfn(1:end-4),'TPWS','pp');
%saveas(h23,fullfile(sdir,ppfn)) 
saveas(h23,fullfile(sdir,ppfn),'png') 

%% Peak Frequency
smsp2 = size(MSP,2);% 2nd element is num fft points
ift = 1:smsp2;
fmsp = ((srate/2)/(smsp2-1))*ift - (srate/2)/(smsp2-1);

[~,im] = max(MSP(:,p.frRange(1):p.frRange(2)),[],2); % maximum between flow-100kHz       
peakFr = fmsp(im + p.frRange(1)-1);

% statistics
mpeakFr = mean(peakFr);
sdpeakFr = std(peakFr);
mepeakFr = median(peakFr);
mopeakFr = mode(peakFr);

% Plot histogram
h24 = figure(24);
nbinsfr = (p.frRange(1):p.frRange(2));
[y,centers] = hist(peakFr,nbinsfr);
bar(centers,y);
xlim([0,srate/2])
title(sprintf('N=%d',length(peakFr)));
xlabel('Peak Frequency (kHz)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpeakFr);
stdlabel = sprintf('Std = %0.2f', sdpeakFr);
melabel = sprintf('Median = %0.2f', mepeakFr);
molabel = sprintf('Mode = %0.2f', mopeakFr);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
% save ici data and figure
pffn = strrep(detfn(1:end-4),'TPWS','peak');
%saveas(h24,fullfile(sdir,pffn))
saveas(h24,fullfile(sdir,pffn),'png')

% %% Excel
% % xls for Danielle
% vecTimes = datevec(MTT(1:end-1));
% icixls = [vecTimes, ici];
% %xlswrite(fn3,icixls);  % write time and click count by bin data to XLS
% % Open Excel, add workbook, change active worksheet,
% % get/put array, save, and close
% % First open an Excel Server
% Excel = actxserver('Excel.Application');
% set(Excel, 'Visible', 1);
% % Insert a new workbook
% Workbooks = Excel.Workbooks;
% Workbook = invoke(Workbooks, 'Add');
% % Make the second sheet active
% Sheets = Excel.ActiveWorkBook.Sheets;
% sheet2 = get(Sheets, 'Item', 1);
% invoke(sheet2, 'Activate');
% % Get a handle to the active sheet
% Activesheet = Excel.Activesheet;
% % Put a MATLAB array into Excel
% %A = [1 2; 3 4];  
% A = icixls;
% la = length(A);
% rstrg = ['A1:G',num2str(la)];
% ActivesheetRange = get(Activesheet,'Range',rstrg);
% set(ActivesheetRange, 'Value', A);
% % Get back a range.  It will be a cell array, 
% % since the cell range can
% % contain different types of data.
% % Range = get(Activesheet, 'Range', 'A1:B2');
% % B = Range.value;
% % % Convert to a double matrix.  The cell array must contain only scalars.
% % B = reshape([B{:}], size(B));
% % Now save the workbook
% xlsfn = strrep(icifn,'.mat','.xls');
% invoke(Workbook, 'SaveAs', xlsfn);
% % To avoid saving the workbook and being prompted to do so,
% % uncomment the following code.
% Workbook.Saved = 1;
% invoke(Workbook, 'Close');
% % Quit Excel
% invoke(Excel, 'Quit');
% % End process
% delete(Excel);

