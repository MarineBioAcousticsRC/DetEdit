function DensityClickGroup(~)
%Density based on click counting and group counting
%JAH 5-10-14, 11/9/2014, 1/2/2016, 2/23/2016, 1/17/2017
% Version makes 4 Figures:
% Fig 1 = click rate by week
% Fig 2 - Group % occur by week
% Fig 3 - Density from click counting
% Fig 4 - Density from group counting
% Define parameters for Density Estimation
% CLICK
% ck = percent false clicks; cvck = CV
% r = click rate #/s; cvr = CV
% w = max range detection; pk = prob detection, cvpk = CV
% GROUP
% ckg = percent false group; cvckg = CV;
% size = group size; cvsize = CV;
% pvg = prob group vocal ; cvpvg = CV;
% wg = max range detection;  pkg =  prob detection; cvpkg = CV;
clear all; % clear data
clf; % clear all figures
nsites = 1; %number of sites to add to the plot
for iplot = 1 : nsites
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
    %     itnum = input('Enter Iteration number (1 2 ...): ','s');
    sp = input('Enter Species: Zc Me BWG Md Ko De Po ','s');
    if (strcmp(sp,'Ko') || strcmp(sp,'k'))
        spe = 'Kogia';
        bindur = 1; % i min bins
        if (strcmp(stn,'GOM'))
            if (strcmp(site,'MC'))
                % density parameters for MC click
                ck = 13.5/100; cvck = 0.11;
                r = 1.41;  cvr = 0.17;
                w = 1.0; pk = 0.0102; cvpk = 0.46/(1.02*sqrt(500)); %210-220
                %BEFORE 08 with 205-210 sl
                % w = 1.0;  pk = 0.0083;  cvpk = 0.42/(0.83*sqrt(500));
                %AFTER 08 with 205-210 sl
                % w = 1.0;  pk = 0.0257;  cvpk = 0.79/(2.57*sqrt(500));
                % density parameters for MC group
                ckg = 1.0/100; cvckg = 0.1;
                size = 2.33; cvsize = 0.15;
                pvg = 0.254; cvpvg = 0.17;
                wg = 1.0; pkg = 0.4254; cvpkg = 4.15/(42.54*sqrt(500));%210-220
                %BEFORE 08 with 205-210 sl
                %  wg = 1.0;  pkg = 0.4819; cvpkg = 4.91/(48.19*sqrt(500));
                %AFTER 08 with 205-210 sl
                %  wg = 1.0;  pkg = 0.5221; cvpkg = 4.55/(52.21*sqrt(500));
            end
            if (strcmp(site,'GC'))
                % density parameters for MC click
                ck = 13.5/100; cvck = 0.11;
                r = 1.41;  cvr = 0.17;
                w = 1.0; pk = 0.0031; cvpk = 0.15/(0.31*sqrt(500));
                %w = 1.0; pk = 0.0033; cvpk = 0.19/(0.33*sqrt(500));%sl205-210
                % density parameters for Mc group
                ckg = 1.0/100; cvckg = 0.1;
                size = 2.33; cvsize = 0.15;
                pvg = 0.254; cvpvg = 0.17;
                wg = 1.0;  pkg = 0.3479; cvpkg = 4.82/(34.79*sqrt(500));
                % wg = 1.0;  pkg = 0.3585; cvpkg = 5.99/(35.85*sqrt(500));%sl205-210
            end
            if (strcmp(site,'DT'))
                % density parameters for MC click
                ck = 13.5/100; cvck = 0.11;
                r = 1.41;  cvr = 0.17;
                w = 1.0; pk = 0.0030; cvpk = 0.16/(0.30*sqrt(500));
                % w = 1.0; pk = 0.0032; cvpk = 0.18/(0.32*sqrt(500));%sl205-210
                % density parameters for Mc group
                ckg = 1.0/100; cvckg = 0.1;
                size = 2.33; cvsize = 0.15;
                pvg = 0.254; cvpvg = 0.17;
                wg = 1.0;  pkg = 0.3461; cvpkg = 4.69/(34.61*sqrt(500));
                % wg = 1.0;  pkg = 0.361; cvpkg = 5.94/(36.1*sqrt(500));%sl205-210
            end
        end
    elseif (strcmp(sp,'Zc') || strcmp(sp,'z'))
        spe = 'Cuviers';
        if (strcmp(stn,'SOCAL'))
            bindur = 1; % 5 min bins
            % density parameters for SOCAL Cuviers click
            ck = 3.0/100; cvck = 0.17;
            r = 0.503; cvr = 0.088;
            w = 4.0;  pk = 0.070;  cvpk = 0.160;
            % density parameters for SOCAL Cuviers group
            ckg = 0.1/100; cvckg = 0.1;
            size = 2.4; cvsize = 0.09;
            pvg = 0.471; cvpvg = 0.09;
            wg = 4.0;  pkg = 0.359; cvpkg = 0.08;
        end
        if (strcmp(stn,'GOM'))
            bindur = 5; % 5 min bins
            if (strcmp(site,'MC'))
                % density parameters for MC Cuviers
                ck = 6.0/100; cvck = 0.04;
                r = 0.493; cvr = 0.088;
                w = 4.0;  pk = 0.070;  cvpk = 0.160;
                % density parameters for MC Cuviers
                ckg = 1.3/100; cvckg = 0.17;
                size = 2.10; cvsize = 0.09;
                pvg = 0.471; cvpvg = 0.09;
                wg = 4.0;  pkg = 0.359; cvpkg = 0.078;
            end
            if (strcmp(site,'GC'))
                % density parameters for GC Cuviers
                ck = 5.6/100; cvck = 0.04;
                r = 0.470; cvr = 0.088;
                w = 4.0;  pk = 0.069; cvpk = 0.158;
                % density parameters for GC Cuviers
                ckg = 0.8/100; cvckg = 0.17;
                size = 1.69; cvsize = 0.10;
                pvg = 0.471; cvpvg = 0.09;
                wg = 4.0;  pkg = 0.360; cvpkg = 0.081;
            end
            if (strcmp(site,'DT'))
                % density parameters for DT Cuviers
                ck = 5.5/100; cvck = 0.04;
                r = 0.457; cvr = 0.087;
                w = 4.0;  pk = .070; cvpk = 0.163;
                % density parameters for DT Cuviers
                ckg = 0.3/100; cvckg = 0.17;
                size = 1.98; cvsize = 0.07;
                pvg = 0.471; cvpvg = 0.09;
                wg = 4.0;  pkg = 0.358; cvpkg = 0.078;
            end
        end
    elseif (strcmp(sp,'Me') || strcmp(sp,'m'))
        spe = 'Gervais';
        if (strcmp(stn,'GOM'))
            bindur = 5; % 5 min bins
            if (strcmp(site,'MC'))
                % density parameters for MC Gervais
                ck = 7.3/100; cvck = 0.04;
                r = 0.492; cvr = 0.169;
                w = 4.0;  pk = 0.043;  cvpk = 0.180;
                % density parameters for MC Gervais
                ckg = 0.7/100; cvckg = 0.17;
                size = 2.18; cvsize = 0.06;
                pvg = 0.254; cvpvg = 0.17
                wg = 4.0;  pkg = 0.281; cvpkg = 0.085;
            end
            if (strcmp(site,'GC'))
                % density parameters for GC Gervais
                ck = 6.4/100; cvck = 0.04;
                r = 0.484; cvr = 0.169;
                w = 4.0;  pk = 0.043; cvpk = 0.178;
                % density parameters for GC Gervais
                ckg = 0.3/100; cvckg = 0.17;
                size = 2.06; cvsize = 0.05;
                pvg = 0.254; cvpvg = 0.17;
                wg = 4.0;  pkg = 0.281; cvpkg = 0.082;
            end
            if (strcmp(site,'DT'))
                % density parameters for DT Gervais
                ck = 4.9/100; cvck = 0.04;
                r = 0.488;  cvr = 0.169;
                w = 4.0;  pk = 0.044;  cvpk = 0.183;
                % density parameters for DT Gervais
                ckg = 0.5/100; cvckg = 0.17;
                size = 2.80; cvsize = 0.08;
                pvg = 0.254; cvpvg = 0.17;
                wg = 4.0;  pkg = 0.278; cvpkg = 0.086;
            end
        end
    elseif (strcmp(sp,'De') || strcmp(sp,'de'))
        spe = 'Delphin';
    elseif (strcmp(sp,'Po') || strcmp(sp,'p'))
        spe = 'Porpoise';
    else
        disp(' Bad Species type')
        return
    end
    % Get Directory with Detections
    disp('Select Directory with Detections');
    sdir = uigetdir('I:\','Select Directory with Detections');
    crfn = [stn,'_',site,'_',spe,'_rate-click.txt'];
    crfigfn = [stn,'_',site,'_',spe,'_rate-click.pdf'];
    fn1t = fullfile(sdir,'\',crfn);    %file with datenum, density
    fn1f = fullfile(sdir,'\',crfigfn);    %file name for pdf figure%
    cdenfn = [stn,'_',site,'_',spe,'_den-click.txt'];
    cdenfigfn = [stn,'_',site,'_',spe,'_den-click.pdf'];
    cdenfigfns = [stn,'_',site,'_',spe,'_den-click-season'];
    fn3t = fullfile(sdir,'\',cdenfn);    %file with datenum, density
    fn3f = fullfile(sdir,'\',cdenfigfn);    %file name for pdf figure%
    fn3fs = fullfile(sdir,'\',cdenfigfns);    %file name for pdf figure%
    grfn = [stn,'_',site,'_',spe,'_rate-group.txt'];
    grfigfn = [stn,'_',site,'_',spe,'_rate-group.pdf'];
    fn2t = fullfile(sdir,'\',grfn);    %file with datenum, density
    fn2f = fullfile(sdir,'\',grfigfn);    %file name for pdf figure%
    gdenfn = [stn,'_',site,'_',spe,'_den-group.txt'];
    gdenfigfn = [stn,'_',site,'_',spe,'_den-group.pdf'];
    gdenfigfns = [stn,'_',site,'_',spe,'_den-group-season'];
    fn4t = fullfile(sdir,'\',gdenfn);    %file with datenum, density
    fn4f = fullfile(sdir,'\',gdenfigfn);    %file name for pdf figure%
    fn4fs = fullfile(sdir,'\',gdenfigfns);    %file name for pdf figure%
    tday = [];
    for n = dpnf:dpnl
        if (site == 'N')% for SOCAL N 31:60
            if ( n ~= 39 && n ~= 42 && n ~= 43 && n ~= 50 && n ~= 55 )
                dpn = num2str(n);
                detpn = [sdir,'\',stn,'_',dpn,'_',site,'\'];
                dayfn = [stn,dpn,site,'_',spe,'_day.txt'];
                fn = fullfile(detpn,dayfn);
                ldata = load(fn,'-ascii');
                tday = [tday;ldata];
            end
        elseif (site == 'H')% for SOCAL H 18:59
            if ( n ~= 19 && n ~= 20 && n ~= 21 && n ~= 22 && n ~= 23 ...
                    && n ~= 24 && n ~= 25 && n ~= 28 && n ~= 30 && n ~= 33 ...
                    && n ~= 39 && n ~= 42 &&n ~= 43 && n ~= 49 ...
                    && n ~= 57);
                dpn = num2str(n);
                detpn = [sdir,'\',stn,'_',dpn,'_',site,'\'];
                dayfn = [stn,dpn,site,'_',spe,'_day.txt'];
                fn = fullfile(detpn,dayfn);
                ldata = load(fn,'-ascii');
                tday = [tday;ldata];
            end
        elseif (site == 'MC')% GOM MC 1:11
            if ( n ~= 8);
                dpn = num2str(n);
                if (n < 10)
                    dpn = ['0',dpn];
                end
                detpn = [sdir,'\'];
                dayfn = [site,dpn,'_',spe,'_day.txt'];
                fn = fullfile(detpn,dayfn);
                ldata = load(fn,'-ascii');
                tday = [tday;ldata];
            end
        elseif (site == 'GC')
            dpn = num2str(n);
            if (n < 10)
                dpn = ['0',dpn];
            end
            detpn = [sdir,'\'];
            dayfn = [site,dpn,'_',spe,'_day.txt'];
            fn = fullfile(detpn,dayfn);
            ldata = load(fn,'-ascii');
            tday = [tday;ldata];
        elseif (site == 'DT')
            dpn = num2str(n);
            if (n < 10)
                dpn = ['0',dpn];
            end
            detpn = [sdir,'\'];
            dayfn = [site,dpn,'_',spe,'_day.txt'];
            fn = fullfile(detpn,dayfn);
            ldata = load(fn,'-ascii');
            tday = [tday;ldata];
        end
    end
    %convert binsDur into sec - add effort as 5th col
    tday(:,5) = bindur*60*tday(:,4); % make 5th column with sec
    % Find repeated tday and add data
    rtday = find(diff(tday(:,1))==0);
    tday(rtday,2:5) = tday(rtday,2:5) + tday(rtday+1,2:5);
    tday(rtday+1,2:5) = tday(rtday,2:5) + tday(rtday+1,2:5);
    [utday,ia,ic] = unique(tday(:,1));
    tday = tday(ia,:);
    if (strcmp(stn,'SOCAL'))
        %         tday(:,1) = tday(:,1) + tj10 + 732; % shift 2009 to 2007
        tday(:,1) = tday(:,1) + tj10 ; % 2009
    elseif (strcmp(stn,'GOM'))
        tday(:,1) = tday(:,1) + tj10;
    end
    %Divide the data into weeks
    % make financial time series fints(dates,data)
    wday = fints(tday);
    mday = tomonthly(wday);
    % EWO = 1 means use Saturday as end of week
    % CumSum returns the cumulative sum of values for each week
    sumweek = toweekly(wday,'CalcMethod','CumSum','BusDays',0,'EOW',1);
    tweekx = getfield(sumweek,'dates');
    tweeky = getfield(sumweek,'series2'); % get # clicks
    tweekyg = getfield(sumweek,'series1'); % get # bins
    stweeky = sum(tweeky);
    stweekyg = sum(tweekyg);
    tweekz = getfield(sumweek,'series4'); % get # sec effort
    tweekzg = getfield(sumweek,'series3'); % get # bins effort
    mtweekz = max(tweekz);
    mtweekzg = max(tweekzg);
    stweekz = sum(tweekz);
    stweekzg = sum(tweekzg);
    % CumSum for month
    summonth = tomonthly(wday,'CalcMethod','CumSum','BusDays',0);
    tmonthx = getfield(summonth,'dates');
    tmonthy = getfield(summonth,'series2'); % get # clicks
    tmonthyg = getfield(summonth,'series1'); % get # bins
    stmonthy = sum(tmonthy);
    stmonthyg = sum(tmonthyg);
    tmonthz = getfield(summonth,'series4'); % get # sec effort
    tmonthzg = getfield(summonth,'series3'); % get # bins effort
    mtmonthz = max(tmonthz);
    mtmonthzg = max(tmonthzg);
    stmonthz = sum(tmonthz);
    stmonthzg = sum(tmonthzg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get each day of week
    fri = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',0);
    sat = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',1);
    sun = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',2);
    mon = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',3);
    tue = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',4);
    wed = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',5);
    thr = toweekly(wday,'CalcMethod','Exact','BusDays',0,'EOW',6);
    tfri = getfield(fri,'series2'); efri = getfield(fri,'series4');
    tsat = getfield(sat,'series2'); esat = getfield(sat,'series4');
    tsun = getfield(sun,'series2'); esun = getfield(sun,'series4');
    tmon = getfield(mon,'series2'); emon = getfield(mon,'series4');
    ttue = getfield(tue,'series2'); etue = getfield(tue,'series4');
    twed = getfield(wed,'series2'); ewed = getfield(wed,'series4');
    tthr = getfield(thr,'series2'); ethr = getfield(thr,'series4');
    tfrig = getfield(fri,'series1'); efrig = getfield(fri,'series3');
    tsatg = getfield(sat,'series1'); esatg = getfield(sat,'series3');
    tsung = getfield(sun,'series1'); esung = getfield(sun,'series3');
    tmong = getfield(mon,'series1'); emong = getfield(mon,'series3');
    ttueg = getfield(tue,'series1'); etueg = getfield(tue,'series3');
    twedg = getfield(wed,'series1'); ewedg = getfield(wed,'series3');
    tthrg = getfield(thr,'series1'); ethrg = getfield(thr,'series3');
    % find day of week for start of time series
    if (wday.dates(1) == sun.dates(1))
        %disp(' Start Sun');
    else
        for i = (length(sun) - 1): -1 : 1
            tsun(i+1) = tsun(i);  esun(i+1) = esun(i);
            tsung(i+1) = tsung(i);  esung(i+1) = esung(i);
        end
        tsun(1) = 0; esun(1) = 0;
        tsung(1) = 0; esung(1) = 0;
        if (wday.dates(1) == mon.dates(1))
            disp(' Start Mon');
        else
            for i = (length(mon) - 1): -1 : 1
                tmon(i+1) = tmon(i); emon(i+1) = emon(i);
                tmong(i+1) = tmong(i); emong(i+1) = emong(i);
            end
            tmon(1) = 0; emon(1) = 0;
            tmong(1) = 0; emong(1) = 0;
            if (wday.dates(1) == tue.dates(1))
                disp(' Start Tue');
            else
                for i = (length(tue) - 1): -1 : 1
                    ttue(i+1) = ttue(i); etue(i+1) = etue(i);
                    ttueg(i+1) = ttueg(i); etueg(i+1) = etueg(i);
                end
                ttue(1) = 0; etue(1) = 0;
                ttueg(1) = 0; etueg(1) = 0;
                if (wday.dates(1) == wed.dates(1))
                    disp(' Start Wed');
                else
                    for i = (length(wed) - 1):  -1 : 1
                        twed(i+1) = twed(i);ewed(i+1) = ewed(i);
                        twedg(i+1) = twedg(i);ewed(i+1) = ewedg(i);
                    end
                    twed(1) = 0; ewed(1) = 0;
                    twedg(1) = 0; ewedg(1) = 0;
                    if (wday.dates(1) == thr.dates(1))
                        disp(' Start Thr');
                    else
                        for i = (length(thr) - 1):  -1 : 1
                            tthr(i+1) = tthr(i); ethr(i+1) = ethr(i);
                            tthrg(i+1) = tthrg(i); ethrg(i+1) = ethrg(i);
                        end
                        tthr(1) = 0; ethr(1) = 0;
                        tthrg(1) = 0; ethrg(1) = 0;
                        if (wday.dates(1) == fri.dates(1))
                            disp(' Start Fri');
                        else
                            for i = (length(fri) - 1):  -1 : 1
                                tfri(i+1) = tfri(i); efri(i+1) = efri(i);
                                tfrig(i+1) = tfrig(i); efrig(i+1) = efrig(i);
                            end
                            tfri(1) = 0; efri(1) = 0;
                            tfrig(1) = 0; efrig(1) = 0;
                            if (wday.dates(1) == sat.dates(1))
                                disp(' Start Sat');
                            end
                        end
                    end
                end
            end
        end
    end
    %make all the days the same length
    mlen = max([length(sun),length(mon),length(tue),length(wed),...
        length(thr),length(fri),length(sat)]);
    if (length(sun) < mlen)
        tsun = [tsun ; 0];  esun = [esun ; 0];
        tsung = [tsung ; 0];  esung = [esung ; 0];
    end
    if (length(mon) < mlen)
        tmon = [tmon ; 0];  emon = [emon ; 0];
        tmong = [tmong ; 0];  emong = [emong ; 0];
    end
    if (length(tue) < mlen)
        ttue = [ttue ; 0];  etue = [etue ; 0];
        ttueg = [ttueg ; 0];  etueg = [etueg ; 0];
    end
    if (length(wed) < mlen)
        twed = [twed ; 0];  ewed = [ewed ; 0];
        twedg = [twedg ; 0];  ewedg = [ewedg ; 0];
    end
    if (length(thr) < mlen)
        tthr = [tthr ; 0];  ethr = [ethr ; 0];
        tthrg = [tthrg ; 0];  ethrg = [ethrg ; 0];
    end
    if (length(fri) < mlen)
        tfri = [tfri ; 0];  efri = [efri ; 0];
        tfrig = [tfrig ; 0];  efrig = [efrig ; 0];
    end
    if (length(sat) < mlen)
        tsat = [tsat ; 0];  esat = [esat ; 0];
        tsatg = [tsatg ; 0];  esatg = [esatg ; 0];
    end
    tdow = [tsun, tmon, ttue, twed, tthr, tfri, tsat];
    edow = [esun, emon, etue, ewed, ethr, efri, esat];
    % calculate weekly variance of Nk/Tk
    te = tdow./edow;
    vte = var(te,0,2);  %any incomplete week is a NaN
    xnan = ~isnan(te);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLICK METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CV for Click Method week
    %cv2nk = vtep'./(tweeky.*tweeky);  %Coeff Var^2 for Nk/Tk
    cv2ck = (cvck * cvck) * ones(length(tweeky),1); % Cv for false click
    cv2pk = (cvpk.*cvpk)* ones(length(tweeky),1); % 0.008 Cv for Probdet function
    cv2r = (cvr.*cvr)* ones(length(tweeky),1); % 0.1 CV for r
    scv2 = cv2ck + cv2pk + cv2r ;  % assume no error Nk/Tk
    % CV for Click Method month
    cv2ckm = (cvck * cvck) * ones(length(tmonthy),1); % Cv for false click
    cv2pkm = (cvpk.*cvpk)* ones(length(tmonthy),1); % 0.008 Cv for Probdet function
    cv2rm = (cvr.*cvr)* ones(length(tmonthy),1); % 0.1 CV for r
    scv2m = cv2ckm + cv2pkm + cv2rm ;  % assume no error Nk/Tk
    %Open Figure for CLICK
    h1 = figure(1);
    %Range of Dates plotted
%     xData = [1,... %2007
%         365+1,... %2008 leap year
%         731+1,... %2009
%         1096+1,... % 2010
%         1461+1,... %2011
%         1826+1,...%2012 leap year
%         2192+1, ...%2013
%         2557+1]; %... %,2014
    %         2922+1, ... %,2015
    %         3287+1, ... %,2016 leap year
    %         3653+1]; % 2017
    %     xData = [731+2-731,... %2009
    %         1096+1-731,... % 2010
    %         1461+1-731,... %2011
    %         1826+1-731+1,...%2012 leap year
    %         2192+1-731+1, ...%2013
    %         2557+1-731, ... %,2014
    %         2922+1-731, ... %,2015
    %         3287+1-731, ... %,2016 leap year
    %         3653+1-731]; % 2017
    xData = [121,244,... %2010
        365+1,365+121,365+244,... %2011
        730+1,730+122,730+245,... %2012 is leap year
        1096+1,1096+121,1096+244];%,... % 2013
%     if (strcmp(stn,'SOCAL'))
        xData = tj10 + xData - 1; % SOCAL
%     elseif (strcmp(stn,'GOM'))
%         xData = tj10 + xData(4:end) -1096 -1; % GOM start 2010
%     end
    mnxD = xData(1);     mxxD = xData(end);
    xTLabel = num2str(xData(1:end-1)');
    %Plot Daily clicks per sec as bar and Weekly as Red Dot %%%%%%%%%%%%
    s(1) = subplot(3,1,iplot); % count
    t3t5 = tday(:,3)./tday(:,5);  % number clicks per sec
    % Detections Plot Scale
    mt3t5 = roundsd(max(1.05*t3t5),2,'ceil');
    yData = 0 : mt3t5/4 : mt3t5;
    yData = roundsd(yData,2,'floor');
    yTLabel = num2str(yData');
    yData(1)= -.05*mt3t5;
    % Determine dates without effort for gray shading
    bpdate = min(xData); epdate = max(xData);
    noeff = []; %dates without effort
    ida = 1;
    cdate = bpdate;
    while (cdate <= epdate)
        if (tday(ida,1) == cdate)
            if (ida < length(tday))
                ida = ida + 1;
            end
        else
            %noeff = [noeff; cdate,mt3t5];
            noeff = [noeff; cdate];
            
        end
        cdate = cdate + 1;
    end
    disp('CLICK METHOD');
    disp([' MeanDaily Clicks/sec=',num2str(mean(t3t5)),...
        ' StDev Clicks/sec= ',num2str(std(t3t5))]);
    disp([' MeanDaily Clicks=',num2str(mean(t3t5*60*60*24)),...
        ' StDev Clicks= ',num2str(std(t3t5*60*60*24))]);
    H1 = bar(tday(:,1),t3t5);
    set(H1,'FaceColor',[1,1,1].*0,'EdgeColor',[1,1,1].*0);
    hold on;
    % for weeks
    tytz = zeros(length(tweekx),1);
    tx = find(tweekz > 1*mtweekz/7);  %Use values with more than 1 day effort
    tytz(tx) = tweeky(tx)./tweekz(tx);  %Weekly average Clicks per sec
    % for months
    tytzm = zeros(length(tmonthx),1);
    txm = find(tmonthz > 1*mtmonthz/30);  %Use values with more than 1 day effort
    tytzm(txm) = tmonthy(txm)./tmonthz(txm);  %Weekly average Clicks per sec
    %
    disp([' MeanWeekly Clicks/sec=',num2str(mean(tytz)),...
        ' StDev Clicks/sec= ',num2str(std(tytz))]);
    disp([' MeanWeekly Clicks=',num2str(mean(tytz*60*60*24*7)),...
        ' StDev Clicks= ',num2str(std(tytz*60*60*24*7))]);
    stytz = stweeky/stweekz; % Entire dataset Clicks per sec
    axis([min(xData),max(xData),min(yData),max(yData)]);
    % Noeffort Bars
    noeff(:,2) = max(yData);
    H2 = bar(noeff(:,1),noeff(:,2));
    set(H2,'FaceColor',[1,1,1]*0.8,'EdgeColor',[1,1,1]*0.8);
    % -3.5 is to shift from Saturday to middle of week
    tweekxp = tweekx - 3.5;
    tmonthxp = tmonthx -15 ; %middle of month
    %plot(tweekxp(tx),tytz(tx),'or','MarkerSize',5,'MarkerFaceColor','r');
    set(s(1),'YTickLabel',yTLabel,'YTick',yData);
    set(s(1),'XTickLabel',xTLabel,'XTick',xData);
    dateFormat = 12;
    datetick('x',dateFormat,'keepticks');
    ylabel('Clicks / sec');
    xticklabel_rotate([],45,[])
    grid on;
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Density Plot for CLICK method
    h3 = figure(3);
    s(1) = subplot(3,1,iplot);
    subplot(s(1));
    % Calculate Density for weekly data
    density = -1 * ones(length(tweekx),1); % when < 1 no effort
    density = [tweekxp,density];
    density(tx,2)=(tytz(tx)*1000*(1-ck))/(pi*w*w*pk*r); %weekly estimate
    sdensity = (stytz*1000*(1-ck))/(pi*w*w*pk*r); %average density
    % Calculate Density Error bar for week
    vden = density(:,2).*density(:,2).*scv2;
    vsden = sdensity * sdensity * scv2(1);
    density = [density, sqrt(vden)];
    sdensityerr = sqrt(vsden);
    % DeSeason
    [deseasTrend,annualChange] = Deseason(stn,site,spe,fn3fs,density);
    figure(3)
    % Density plot for month
    densitym = -1 * ones(length(tmonthx),1); % when < 1 no effort
    densitym = [tmonthxp,densitym];
    densitym(txm,2)=(tytzm(txm)*1000*(1-ck))/(pi*w*w*pk*r); %weekly estimate
    sdensitym = (stytz*1000*(1-ck))/(pi*w*w*pk*r); %average density
    % Calculate Density Error bar for week
    vdenm = densitym(:,2).*densitym(:,2).*scv2m;
    vsdenm = sdensitym * sdensitym * scv2m(1);
    densitym = [densitym, sqrt(vdenm)];
    sdensityerrm = sqrt(vsdenm);
    % Determine Plot Y scale
    mdensity = roundsd(max(1.15*density(:,2)),2,'ceil');
    yData = 0 : mdensity/4 : mdensity;
    yTLabel = num2str(yData');
    yData(1)= -.05*mdensity;
    % Noeffort Bars
    noeff(:,2) = max(yData);
    H3 = bar(noeff(:,1),noeff(:,2));
    hold on;
    set(H3,'FaceColor',[1,1,1]*0.8,'EdgeColor',[1,1,1]*0.8);
    % Density dots
    H4 = errorbar(tweekxp(tx),density(tx,2),density(tx,3),'ok',...
        'MarkerSize',4,'MarkerFaceColor','w');
    errorbar_tick(H4,0,'UNITS');
    % H5 = errorbar(tmonthxp(txm),densitym(txm,2),density(txm,3),'or',...
    %     'MarkerSize',6,'MarkerFaceColor','w');
    % errorbar_tick(H5,0,'UNITS');
    %% Plot Regression Line
    %     plot(tweekxp,deseasTrend,'--','LineWidth',2.0);
    %     title(['Annual Trend = ',num2str(annualChange)]);
    set(s(1),'XTick',xData);
    set(s(1),'XLim',[mnxD, mxxD]);
    set(s(1),'YLim',[min(yData), max(yData)]);
    set(s(1),'YTickLabel',yTLabel,'YTick',yData);
    dateFormat = 12;
    datetick('x',dateFormat,'keepticks')
    %xlabel('Date');
    ylabel('Density #/1000 km^2');
    grid on;
    hold off;
    xticklabel_rotate([],45,[])
    %subplot(3,1,iplot)
    % slope = 365*round(10000*b(2))/10000;
    % slopepm = 365*round(10000*(bint(2,2)-bint(2,1)))/(2*10000);
    % tstr = {[stn,' ',spe,' Rate= ',num2str(stytz),' C/s',...
    %     '  D= ',num2str(sdensity),'+- ',num2str(sdensityerr)]};
    % %     ' Sl=',num2str(slope),'+- ',num2str(slopepm)]};
    %title(tstr);
    tstr = {[stn,' ',site,' ',spe,' Rate= ',num2str(stytz),' C/s',...
        '  D= ',num2str(sdensity),'+- ',num2str(sdensityerr)]};
    title(tstr);
    %     save(fn1,'density','-ascii');
    %     print(h,'-painters','-dpdf','-r600',fn2);
    % set(s2,'XTick',xData);
    % set(s2,'XLim',[mnxD, mxxD]);
    % set(s2,'YTickLabel',yTLabel,'YTick',yData);
    hold off;
    disp([' Sum Rate= ',num2str(stytz),' Click/s',...
        '  Density= ',num2str(sdensity),' Err =',num2str(sdensityerr)]);
    %     ' Slope=',num2str(slope),...
    %     ' -Slope=',num2str(365*bint(2,1)),...
    %       ' +Slope=',num2str(365*bint(2,2))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GROUP METHOD %%%%%%%%%%%%%%%%%%%%%
    % CV for Group Method
    cv2ckg = (cvckg * cvckg) * ones(length(tweekyg),1); % 0.025 Cv for false click
    cv2pvg = (cvpvg.*cvpvg) * ones(length(tweekyg),1); % Cv for Prob group vocal
    cv2pkg = (cvpkg.*cvpkg)* ones(length(tweekyg),1); % .004 Cv for Probdet function
    cv2sizeg = (cvsize*cvsize)* ones(length(tweekyg),1); % CV for group size
    scv2g = cv2sizeg + cv2pvg + cv2pkg + cv2ckg;
    %Open Figure
    h2 = figure(2);
    % xDATA same as for click
    % Dates without effort for gray shading - same as click
    %Plot Weekly #of bins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s(1) = subplot(3,1,iplot);
    t2t4 = tday(:,2)./tday(:,4);
    % Detections Plot Scale
    mt2t4 = roundsd(max(105*t2t4),2,'ceil'); %percentage
    yData = 0 : mt2t4/4 : mt2t4;
    yData = roundsd(yData,2,'floor');
    yTLabel = num2str(yData');
    yData(1)= -.05*mt2t4;
    disp('GROUP METHOD');
    disp([' MeanDaily frac Groups=',num2str(mean(t2t4)),...
        ' StDev frac Groups= ',num2str(std(t2t4))]);
    disp([' MeanDaily Groups=',num2str(mean(t2t4*12*24)),...
        ' StDev Groups= ',num2str(std(t2t4*12*24))]);
    H1 = bar(tday(:,1),100*t2t4); %percentage
    set(H1,'FaceColor',[1,1,1]*0,'EdgeColor',[1,1,1]*0);
    hold on;
    tytzg = zeros(length(tweekx),1);
    tx = find(tweekz > 1*mtweekz/7);  %Use values with more than 1 day effort
    tytzg(tx) = tweekyg(tx)./tweekzg(tx);  %Tk/Nk for weekly average
    disp([' MeanWeekly frac Groups=',num2str(mean(tytzg(tx))),...
        ' StDev frac Groups= ',num2str(std(tytzg(tx)))]);
    disp([' MeanWeekly Groups=',num2str(mean(tytzg(tx)*12*24*7)),...
        ' StDev Groups= ',num2str(std(tytzg(tx)*12*24*7))]);
    stytzg = stweekyg/stweekzg; % Tk/Nk for entire dataset
    % -3.5 is to shift from Saturday to middle of week
    axis([min(xData),max(xData),min(yData),max(yData)]);
    % no effort gray bars
    noeff(:,2) = max(yData);
    H2 = bar(noeff(:,1),noeff(:,2));
    set(H2,'FaceColor',[1,1,1]*0.8,'EdgeColor',[1,1,1]*0.8);
    tweekxp = tweekx - 3.5;
    %plot(tweekxp(tx),tytz(tx),'or','MarkerSize',5,'MarkerFaceColor','r');
    set(s(1),'YTickLabel',yTLabel,'YTick',yData);
    set(s(1),'XTick',xData);
    dateFormat = 12;
    datetick('x',dateFormat,'keepticks');
    binlab = num2str(bindur);
    ylabel(['Daily %',binlab,' min bins']);
    xticklabel_rotate([],45,[])
    grid on;
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Density Plot
    h4 = figure(4);
    s(1) = subplot(3,1,iplot);
    subplot(s(1));
    % Calculate Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    densityg = -1 * ones(length(tweekx),1); % when < 1 no effort
    densityg = [tweekxp,densityg];
    densityg(tx,2)=(tytzg(tx)*1000*size*(1-ckg))/(pi*wg*wg*pkg*pvg); %weekly estimate
    sdensityg = (stytzg*1000*size*(1-ckg))/(pi*wg*wg*pkg*pvg); %average density
    % Calculate Density Error bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vdeng = densityg(:,2).*densityg(:,2).*scv2g;
    vsdeng = sdensityg * sdensityg * scv2g(1);
    densityg = [densityg, sqrt(vdeng)];
    sdensityerrg = sqrt(vsdeng);
    % DeSeason
    [deseasTrend,annualChange] = Deseason(stn,site,spe,fn4fs,densityg);
    figure(4)
    % Determine Plot Y scale
    mdensityg = roundsd(max(1.2*densityg(:,2)),2,'ceil');
    yData = 0 : mdensityg/4 : mdensityg;
    yData = roundsd(yData,2,'floor');
    yTLabel = num2str(yData');
    yData(1)= -.05*mdensityg;
    %No effort Bars
    noeff(:,2) = max(yData);
    H3 = bar(noeff(:,1),noeff(:,2));
    set(H3,'FaceColor',[1,1,1]*0.8,'EdgeColor',[1,1,1]*0.8);
    hold on;
    % Density dots
    H4 = errorbar(tweekxp(tx),densityg(tx,2),densityg(tx,3),'ok',...
        'MarkerSize',4,'MarkerFaceColor','w');
    errorbar_tick(H4,0,'UNITS')
    %%Plot Regression Line
    %     plot(tweekxp,deseasTrend,'--','LineWidth',2.0);
    %     title(['Annual Trend = ',num2str(annualChange)]);
    set(s(1),'XTick',xData);
    set(s(1),'XLim',[mnxD, mxxD]);
    set(s(1),'YLim',[min(yData), max(yData)]);
    set(s(1),'YTickLabel',yTLabel,'YTick',yData);
    dateFormat = 12;
    datetick('x',dateFormat,'keepticks')
    % xlabel('Date');
    ylabel('Density #/1000 km^2');
    grid on;
    xticklabel_rotate([],45,[])
    hold off;
    %subplot(3,1,iplot)
    % slope = 365*round(10000*b(2))/10000;
    % slopepm = 365*round(10000*(bint(2,2)-bint(2,1)))/(2*10000);
    %     tstr = {[stn,' ',spe,' Rate= ',num2str(stytzg),' Bins',...
    %         '  D= ',num2str(sdensityg),'+- ',num2str(sdensityerrg)]};
    %     %    % ' Sl=',num2str(slope),'+- ',num2str(slopepm)]};
    %     title(tstr);
    tstr = {[stn,' ',site,' ',spe,' Rate= ',num2str(stytzg),' % bins',...
        '  D= ',num2str(sdensityg),'+- ',num2str(sdensityerrg)]};
    title(tstr);
    %xticklabel_rotate
    %     save(fn3,'densityg','-ascii');
    %     print(h,'-painters','-dpdf','-r600',fn4);
    % set(s2,'XTick',xData);
    % set(s2,'XLim',[mnxD, mxxD]);
    % set(s2,'YTickLabel',yTLabel,'YTick',yData);
    %hold off;
    disp([' Sum Rate= ',num2str(stytzg),' Fraction Bins',...
        '  Density= ',num2str(sdensityg),' Err =',num2str(sdensityerrg)]);
    % ' Slope=',num2str(slope),...
    %' -Slope=',num2str(365*bint(2,1)),...
    %' +Slope=',num2str(365*bint(2,2))]);
end
%save data and figure % Save data and figure%save data and figure
st3t5 = [tday(:,1),t3t5];
save(fn1t,'st3t5','-ascii');
print(h1,'-painters','-dpdf','-r600',fn1f);
save(fn3t,'density','-ascii');
print(h3,'-painters','-dpdf','-r600',fn3f);
st2t4 = [tday(:,1),t2t4];
save(fn2t,'st2t4','-ascii');
print(h2,'-painters','-dpdf','-r600',fn2f);
save(fn4t,'densityg','-ascii');
print(h4,'-painters','-dpdf','-r600',fn4f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                   Subroutines                           %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=roundsd(x,n,method)
%ROUNDSD Round with fixed significant digits
%	ROUNDSD(X,N) rounds the elements of X towards the nearest number with
%	N significant digits.
%
%	ROUNDSD(X,N,METHOD) uses following methods for rounding:
%		'round' - nearest (default)
%		'floor' - towards minus infinity
%		'ceil'  - towards infinity
%		'fix'   - towards zero
%
%	Examples:
%		roundsd(0.012345,3) returns 0.0123
%		roundsd(12345,2) returns 12000
%		roundsd(12.345,4,'ceil') returns 12.35
% --- but to avoid numerical noise, we must treat separately positive and
% negative exponents, because:
% 3.55/0.1 - 35.5 is -7.105427357601e-15
% 	3.55*10 - 35.5 is 0
e = floor(log10(abs(x)) - n + 1);
og = 10.^abs(e);
y = feval(method,x./og).*og;
k = find(e<0);
if ~isempty(k)
    y(k) = feval(method,x(k).*og(k))./og(k);
end
y(x==0) = 0;

end