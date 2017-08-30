function [ ] = Calicippfunc(MTT,MPP,filePrefix,sp,sdir,detfn,p)
% Calculate the ICI and the PP 
close all

%% Inter-Click Interval
ici = diff(MTT)*24*60*60;

%define threshold for ici
if isempty(p.iciMin)
    p.iciMin = min(ici);
end
if isempty(p.iciMax)
    p.iciMax = max(ici);
end
iciSel = ici(ici > p.iciMin & ici < p.iciMax);

% plot ici histogram
h22 = figure(22);
hist(iciSel,100);
miciSel = mean(iciSel);
sdiciSel = std(iciSel);
title(sprintf('\\mu=%0.2f \\sigma=%0.2f N=%d',miciSel,sdiciSel,length(MTT)))
xlabel(sprintf('Inter-Pulse Interval (s) - thr: %0.2f-%0.2f',p.iciMin,p.iciMax))

% save ici data and figure
icifn = strrep(detfn,'TPWS','ici');
save(fullfile(sdir,icifn),'iciSel')
saveas(h22,fullfile(sdir,icifn),'m')  % save as matlab fig

%% Peak-to-peak
% Plot peak-peak data
h23 = figure(23);
hist(MPP,100);
mpp = mean(MPP);
sdpp = std(MPP);
title(sprintf([filePrefix,' ',sp,' Mean= ',num2str(mpp),...
    'dB  StDev= ',num2str(sdpp),' Number= ',num2str(length(MPP))]))
title(sprintf('\\mu=%0.2f \\sigma=%0.2f N=%d',mpp,sdpp,length(MPP)))
xlabel('Peak-Peak Amplitude')

% save pp figure
ppfn = strrep(detfn,'TPWS','pp');
saveas(h23,fullfile(sdir,ppfn),'m') % save as matlab fig

% xls for Danielle
icixls = datevec(MTT(1:end-1));
icixls = [icixls, ici];
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
la = length(A);
rstrg = ['A1:G',num2str(la)];
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

