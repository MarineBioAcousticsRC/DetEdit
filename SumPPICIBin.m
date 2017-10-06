clear all; % clear data
clf; % clear all figures
fignum = 200;
spe = 'Kogia';
icimin = 0.05; icimax = .3;
stn = input('Enter Project Name (SOCAL or GOM): ','s'); % site name
if (strcmp(stn,'SOCAL'))
    tj10 = datenum([2009 1 1 0 0 0]);   %Jan 2009 is time = 0
elseif (strcmp(stn,'GOM'))
    tj10 = datenum([2010 1 1 0 0 0]);   %Jan 2010 is time = 0
end
site = input('Enter site (N, M, H): ','s');
dpnfs = input('Enter First Deployment number (01 02 ...): ','s');
dpnf = str2double(dpnfs);
dpnls = input('Enter Last Deployment number (01 02 ...): ','s');
dpnl = str2double(dpnls);
disp('Select Directory with Detections');
sdir = uigetdir('I:\','Select Directory with Detections');
pp = []; ici = []; bin =[]; diel =[]; td = [];
for n = dpnf:dpnl
    if (site == 'MC')% GOM MC 1:11
        if ( n ~= 8);
            dpn = num2str(n);
            if (n < 10)
                dpn = ['0',dpn];
            end
        end
    elseif (site == 'GC')
        dpn = num2str(n);
        if (n < 10)
            dpn = ['0',dpn];
        end
    elseif (site == 'DT')
        dpn = num2str(n);
        if (n < 10)
            dpn = ['0',dpn];
        end
    end
    detpn = [sdir,'\'];
    % read Test Detection data
    ppfn = [site,dpn,'_',spe,'_TD4.mat'];
    fn0 = fullfile(detpn,ppfn);
    load(fn0);
    itd = find(zTD(:,1) > 0 & zTD(:,2) > -0.5);
    ntest = sum(zTD(itd,:));
    disp([site,dpn,' Num False = ',num2str(ntest(2)),...
        '  Num Tested = ',num2str(ntest(1)),...
        ' False Percent = ',num2str(100*ntest(2)/ntest(1))]);
    td = [td;zTD];
    % read pp data
    ppfn = [site,dpn,'_',spe,'_pp.mat'];
    fn1 = fullfile(detpn,ppfn);
    ldata = load(fn1);
    x = ldata.ppnf;
    if isrow(x)
        pp = [pp; x'];
    else
        pp = [pp; x];
    end
    % read ici data
    icifn = [site,dpn,'_',spe,'_ici2.xls'];
    fn2 = fullfile(detpn,icifn);
    icidata = xlsread(fn2);
    ici = [ici;icidata(:,7)];
    
    % read Max PP by bin data
    binfn = [site,dpn,'_',spe,'_bin.xls'];
    fn3 = fullfile(detpn,binfn);
    bindata = xlsread(fn3);
    bin = [bin;bindata(:,8)];
    diel = [diel;bindata(:,4)+(bindata(:,5)/60)];
end
%
%Calculate False positive rate
id = find(td(:,1) > 0 & td(:,2) > -0.5);
false = sum(td(id,:));
disp([' Num False = ',num2str(false(2)),...
    '  Num Tested = ',num2str(false(1)),...
    ' False Percent = ',num2str(100*false(2)/false(1))]);
%plot ici
iciIdx = find(ici > icimin & ici < icimax);
figure(fignum); fignum = fignum +1;
hist(ici(iciIdx),200);
%set(gca,'EdgeColor','k','FaceColor','w')
icif = ici(iciIdx);
micif = mean(icif);
sdicif = std(icif);
icistr=[spe,' ',site,' ','Mean= ',num2str(micif),...
    ' StDev= ',num2str(sdicif),' Number= ',num2str(length(iciIdx))];
title(icistr)
xlabel('Inter-Pulse Interval (s)');
icifn = [site,'_',spe,'_ici.pdf'];
fn1 = fullfile(detpn,icifn);
saveas(gcf,fn1,'pdf')

%plot pp click data
figure(fignum); fignum = fignum +1; %linear histogram
%hist(pp,100);
center = 115.5:1:160.5;
lpp = length(pp);
mpp = mean(pp);
sdpp = std(pp);
%[n, center] = hist(pp,100);
[nhist] = hist(pp,center);
h2 = bar(center, nhist, 'barwidth', 1, 'basevalue', 1);
set(h2,'EdgeColor','k','FaceColor','w')
xlim([min(center),max(center)])
ylim([0,1.05*max(nhist)])
ylabel(gca, 'Number of detections','FontSize',16)
xlabel(gca,'Peak-Peak Amplitude (dB)','FontSize',16);
ppstr=[spe,' ',site,' ','Mean= ',num2str(mpp),...
    'dB  StDev= ',num2str(sdpp),' Number= ',num2str(length(pp))];
title(ppstr)

%plot diel pattern
figure(fignum); fignum = fignum +1;
sutc = 5;% shift to local midnight
dm = mod(diel - sutc,24); 
dcenter = 0.5:1:23.5;
[dhist] = hist(dm,dcenter);
hd = bar(dcenter, dhist, 'barwidth', 1, 'basevalue', 1);
set(hd,'EdgeColor','k','FaceColor','w')
xlim([0,24]);
ylim([0,1.05*max(dhist)])
ylabel(gca, 'Number of Bins','FontSize',16)
xlabel(gca,['Time of Day (UTC-',num2str(sutc),')'],'FontSize',16);
dstr=[spe,' ',site,' Number Bins= ',num2str(length(diel))];
title(dstr);
dfn = [site,'_',spe,'_diel.pdf'];
fnd = fullfile(detpn,dfn);
saveas(gcf,fnd,'pdf')

figure(fignum); fignum = fignum +1; %percent log histogram PP click
nper = nhist*100/lpp; % percentage
h3 = bar(center,nper, 'barwidth', 1, 'basevalue', 1);
set(h3,'EdgeColor','k','FaceColor','w')
%plot(center,n*100/lpp);
set(gca,'YScale','log')
xlim([min(center),max(center)])
ylim([.1,50])
set(gca,'FontSize',12)
title(ppstr)
ylabel(gca, 'Percent of detections','FontSize',16)
xlabel(gca,'Peak-Peak Amplitude (dB)','FontSize',16);
pplog = [site,'_',spe,'_pplog.mat'];
fnlog = fullfile(detpn,pplog);
save(fnlog,'center','nper')
ppfn = [site,'_',spe,'_pp.pdf'];
fn3 = fullfile(detpn,ppfn);
saveas(gcf,fn3,'pdf')

%plot pp bin data
figure(fignum); fignum = fignum +1; %linear histogram
%hist(pp,100);
lbin = length(bin);
mbin = mean(bin);
sdbin = std(bin);
%[n, center] = hist(pp,100);
[nbhist] = hist(bin,center);
h2 = bar(center, nbhist, 'barwidth', 1, 'basevalue', 1);
set(h2,'EdgeColor','k','FaceColor','w')
xlim([min(center),max(center)])
ylim([0,1.05*max(nbhist)])
ylabel(gca, 'Number of bins','FontSize',16)
xlabel(gca,'Peak-Peak Amplitude (dB)','FontSize',16);
binstr=[spe,' ',site,' ','Mean= ',num2str(mbin),...
    'dB  StDev= ',num2str(sdbin),' Number= ',num2str(length(bin))];
title(binstr)

figure(fignum); fignum = fignum +1; %percent log histogram PP click
nbper = nbhist*100/lbin; % percentage
h3 = bar(center,nbper, 'barwidth', 1, 'basevalue', 1);
set(h3,'EdgeColor','k','FaceColor','w')
%plot(center,n*100/lpp);
set(gca,'YScale','log')
xlim([min(center),max(center)])
ylim([.1,50])
set(gca,'FontSize',12)
title(binstr)
ylabel(gca, 'Percent of bins','FontSize',16)
xlabel(gca,'Peak-Peak Amplitude (dB)','FontSize',16);
binlog = [site,'_',spe,'_binlog.mat'];
fnblog = fullfile(detpn,binlog);
save(fnblog,'center','nbper')
binfn = [site,'_',spe,'_bin.pdf'];
fn5 = fullfile(detpn,ppfn);
saveas(gcf,fn5,'pdf')
