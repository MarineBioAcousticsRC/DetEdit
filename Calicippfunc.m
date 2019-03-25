function [ ] = Calicippfunc(MTT,MPP,MSP,outDir,detfn,p)
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
figure(1); set(1,'name','Inter-Click Interval')
h1 = gca;
centerIci = (p.iciRange(1):p.iciRange(2));
[nhist] = histc(iciSel,centerIci);
bar(h1,centerIci,nhist, 'barwidth', 1, 'basevalue', 1);
xlim(h1,p.iciRange);
ylim(h1,[1 inf]);
title(h1,sprintf('N=%d',length(MTT)))
xlabel(h1,'Inter-Pulse Interval (ms)')
ylabel(h1,'Counts')
% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', miciSel);
stdlabel = sprintf('Std = %0.2f', sdiciSel);
melabel = sprintf('Median = %0.2f', meiciSel);
molabel = sprintf('Mode = %0.2f', moiciSel);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% save ici data and figure
icifn = strrep(detfn(1:end-4),'TPWS','ici');
% saveas(h1,fullfile(sdir,icifn))
saveas(h1,fullfile(outDir,icifn),'png')

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
figure(2); set(2,'name','Received Levels')
h2 = gca;
center = p.threshRL:1:p.p1Hi;
[nhist] = histc(MPP,center);
bar(h2,center, nhist, 'barwidth', 1, 'basevalue', 1)
xlim(h2,[p.threshRL-0.5, p.p1Hi+0.5])
ylim(h2,[1 inf]);
title(h2,sprintf('N=%d',length(MPP)))
xlabel(h2,'Peak-Peak Amplitude (dB)')

% create labels and textbox
mnlabel = sprintf('Mean = %0.2f', mpp);
stdlabel = sprintf('Std = %0.2f', sdpp);
melabel = sprintf('Median = %0.2f', mepp);
molabel = sprintf('Mode = %0.2f', mopp);
annotation('textbox',[0.58 0.75 0.1 0.1],'String',{mnlabel,stdlabel,...
    melabel,molabel});

% Save plot
ppfn = strrep(detfn(1:end-4),'TPWS','pp');
%saveas(h2,fullfile(sdir,ppfn)) 
saveas(h2,fullfile(outDir,ppfn),'png') 

%% Peak Frequency
smsp2 = size(MSP,2);% 2nd element is num fft points
ift = 1:smsp2;
fmsp = ((p.sampleRate/2)/(smsp2-1))*ift - (p.sampleRate/2)/(smsp2-1);

[~,im] = max(MSP(:,p.frRange(1):p.frRange(2)),[],2); % maximum between flow-100kHz       
peakFr = fmsp(im + p.frRange(1)-1);

% statistics
mpeakFr = mean(peakFr);
sdpeakFr = std(peakFr);
mepeakFr = median(peakFr);
mopeakFr = mode(peakFr);

% Plot histogram
figure(3); set(3,'name','Peak Frequency')
h3 = gca;
nbinsfr = (p.frRange(1):p.frRange(2));
[y,centers] = hist(peakFr,nbinsfr);
bar(h3,centers,y);
xlim(h3,[0,p.sampleRate/2])
title(h3,sprintf('N=%d',length(peakFr)));
xlabel(h3,'Peak Frequency (kHz)')

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
saveas(h3,fullfile(p.sdir,pffn))


