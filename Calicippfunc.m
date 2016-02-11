function [ ] = Calicippfunc(MTT,MPP,stn,dpn,spe,pn1,icimin,icimax)
% Calculate the ICI and the PP 

icifn = [stn,dpn,'_',spe,'_ici'];
icifnfig = [stn,dpn,'_',spe,'_ici2.pdf'];
xlsfn = [stn,dpn,'_',spe,'_ici2.xls'];
ppfn = [stn,dpn,'_',spe,'_pp'];
ppfnfig = [stn,dpn,'_',spe,'_pp2.pdf'];
ppxlsfn = [stn,dpn,'_',spe,'_pp2.xls'];
% fn = fullfile(detpn,detfn);
fn1 = fullfile(pn1,'\',icifn);
fn2 = fullfile(pn1,'\',icifnfig);
fn3 = fullfile(pn1,'\',xlsfn);
fn4 = fullfile(pn1,'\',ppfn);
fn5 = fullfile(pn1,'\',ppfnfig);
%
ctnf = MTT;
ppnf = MPP;
% end
ici = diff(ctnf)*24*60*60;
%define threshold for short ici to be 50 ms and long to be 1000 ms
iciIdx = find(ici > icimin & ici < icimax);
figure(22)
hist(ici(iciIdx),100);
icif = ici(iciIdx);
micif = mean(icif);
sdicif = std(icif);
icistr=[stn,' ',dpn,' ',spe,' ','Mean= ',num2str(micif),...
    ' StDev= ',num2str(sdicif),' Number= ',num2str(length(ctnf))];
title(icistr)
xlabel('Inter-Pulse Interval (s)')
%save(fn1,'icif','-ascii')
save(fn1,'icif')
saveas(gcf,fn2,'pdf')
% Save peak-peak data
figure(23)
hist(ppnf,100);
mpp = mean(ppnf);
sdpp = std(ppnf);
ppstr=[stn,' ',dpn,' ',spe,' ','Mean= ',num2str(mpp),...
    'dB  StDev= ',num2str(sdpp),' Number= ',num2str(length(ppnf))];
title(ppstr)
xlabel('Peak-Peak Amplitude')
save(fn4,'ppnf')
saveas(gcf,fn5,'pdf')
% xls for Danielle
icixls = datevec(ctnf(1:end-1));
icixls = [icixls, ici'];
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

