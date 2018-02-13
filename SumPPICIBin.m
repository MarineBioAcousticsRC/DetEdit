function SumPPICIBin(varargin)
fignum = 200;

% get user input and set up file names
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
        case 'countType'
            countType = varargin{n+1}; n=n+2;
        case 'effort'
            effort = varargin{n+1}; n=n+2;
        case 'refTime'
            refStartTime = varargin{n+1}; n=n+2;
        otherwise
            error('Bad optional argument: "%s"', varargin{n});
    end
end

p = sp_setting_defaults('sp',sp,'analysis','SumPPICIBin');

PPall = []; TTall = [];
ICIall = [];
bdAll = [];
bin =[]; diel =[]; tdAll = [];
for idsk = 1 : length(concatFiles)
    % Load file
    sprintf('Loading file %s',fullfile(sdir,concatFiles{idsk}))
    load(fullfile(sdir,concatFiles{idsk}));
    
    % Group times and peak-to-peak
    PPall = [PPall; MPP];
    TTall = [TTall; MTT];
    
    if countType == 'C' || countType =='B'
        % Inter-Click Interval
        ici = diff(MTT)*24*60*60*1000; % in ms
        ICIall = [ICIall;ici];
    end
end

% exclude times outside of effort
outeff = find(TTall >= effort.Start & TTall <= effort.End);

    if countType == 'G' || countType =='B'
        % convert times to bin vector times
        vTT = datevec(TTall);
        tbin = datetime([vTT(:,1:4), floor(vTT(:,5)/p.binDur)*p.binDur, ...
            zeros(length(vTT),1)]);
        
        % create table and get click counts and max pp per bin
        TTall = MPP; PPall = MPP;
        data = timetable(tbin,TTall,PPall);
        binData = varfun(@max,data,'GroupingVariable','tbin','InputVariable','PPall');
        binData.Properties.VariableNames{'GroupCount'} = 'ClickCount';
        binData.Properties.VariableNames{'max_ppAll'} = 'maxPP';
        
        % group data by days
        dayData = retime(binData,'daily','count');
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create plots, binlog and pplog files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if countType == 'C' || countType =='B'
    %% Plot ici
    iciIdx = find(ICIall > p.iciRange(1) & ICIall < p.iciRange(2));
    % statistics
    miciSel = mean(ICIall(iciIdx));
    sdiciSel = std(ICIall(iciIdx));
    moiciSel = mode(ICIall(iciIdx));
    meiciSel = median(ICIall(iciIdx));
    
    h1 = figure(fignum); fignum = fignum +1;
    nbinsici = (p.iciRange(1):p.iciRange(2));
    [y,centers] = hist(ICIall(iciIdx),nbinsici);
    bar(centers,y);
    xlim(p.iciRange);
    title(sprintf('%s N=%d',filePrefix,length(ICIall)), 'Interpreter', 'none')
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
    icifn = [filePrefix,'_',p.speName,'_ici'];
    saveas(h1,fullfile(sdir,icifn),'png')
    
    %% Plot pp click data
    % statistics
    mpp = mean(PPall);
    sdpp = std(PPall);
    mopp = mode(PPall);
    mepp = median(PPall);
    
    % Plot histogram
    h2 = figure(fignum); fignum = fignum +1; %linear histogram
    center = p.threshRL:1:p.p1Hi;
    [nhist] = histc(PPall,center);
    bar(center, nhist, 'barwidth', 1, 'basevalue', 1)
    title(sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
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
    ppfn = [filePrefix,'_',p.speName,'_pp'];
    saveas(h2,fullfile(sdir,ppfn),'png')
    
    %% Percent log histogram PP click
    h3 = figure(fignum); fignum = fignum +1; %percent log histogram PP click
    nper = nhist*100/length(PPall); % percentage
    bar(center,nper, 'barwidth', 1, 'basevalue', 1,'EdgeColor','k','FaceColor','w');
    set(gca,'YScale','log')
    xlim([min(center)-0.5,max(center)+0.5])
    %ylim([.1,50])
    set(gca,'FontSize',12)
    title(sprintf('%s N=%d',filePrefix,length(PPall)), 'Interpreter', 'none')
    ylabel('Percent of detections','FontSize',16)
    xlabel('Peak-Peak Amplitude (dB)','FontSize',16);
    
    % Save plot and pplog values in a mat file
    pplog = [filePrefix,'_',p.speName,'_pplog'];
    save(fullfile(sdir,pplog),'center','nper')
    saveas(h3,fullfile(sdir,pplog),'png')   
end

if countType == 'G' || countType =='B' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to modify
    
    
    
    
    
    
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
    dstr=[p.speName,' ',site,' Number Bins= ',num2str(length(diel))];
    title(dstr);
    dfn = [site,'_',p.speName,'_diel.pdf'];
    fnd = fullfile(detpn,dfn);
    saveas(gcf,fnd,'pdf')
    
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
    binstr=[p.speName,' ',site,' ','Mean= ',num2str(mbin),...
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
    binlog = [site,'_',p.speName,'_binlog.mat'];
    fnblog = fullfile(detpn,binlog);
    save(fnblog,'center','nbper')
    binfn = [site,'_',p.speName,'_bin.pdf'];
    fn5 = fullfile(detpn,ppfn);
    saveas(gcf,fn5,'pdf')
end

