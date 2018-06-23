% brushDet_171002.m
%
% brush Detection Times via tracks in TDOAs, AZ, EL plots
%
% 171002 smw
clear
global hf
%clearvars
% %
% fstr = 'HAT_B_01_02_C4_170510_1216';
% ifname = [fstr,'_DetCL1.mat'];
% ofname = [fstr,'_DetCL2.mat'];
% pname = 'F:\HAT_B_01_array\matfiles\BW\';
ifname = 'Fig201ses1.fig';
tpws = 'HAT04A_diska_Delphin_TPWS1.mat';
pname = '/Users/jah/DCLDE/HAT_A/';
ffile = fullfile(pname,ifname);
tfile = fullfile(pname,tpws);
% sfile = fullfile(pname,ofname);

DI = {};
%load(dfile)  % Det structure: Time, Ang, RL, TDOA, TDOA2
%
load(tfile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc = 'a';

if isempty(DI) % input not previously brushed
%     Ia = 1:length(Det.Time);    %  original indices of detections
    Ia = 1:length(MTT);    %  original indices of detections JAH
    Ib = 0;                     % 'loser'/brushed away indices
    Ic = Ia;                    % unbrushed 'keeper' indices
else          % input previously brushed
    Ia = DI.C;
    Ib = DI.B;
    Ic = Ia;
end
hf = openfig(ffile);
while ~strcmp(cc,'q')
    %
%     hf = figure(1000);
%     clf
%     fs = 10;
%     ms = 4;
%     
%     subplot(2,1,1)
%     plot(Det.Time(Ic),Det.Ang(Ic,1),'b.','MarkerSize',ms);
%     
%     title(dfile,'Interpreter','none','FontSize',fs)
%     h=gca;
%     set(h,'XMinorTick','on')
%     set(h,'FontSize',fs);
%     ylabel('Azimuth [deg]')
%     grid minor
%     v=axis;
% %     axis([v(1) v(2) 0 360])
%     datetick('x',13)
%     set(h,'XTickLabel',[])
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(2,1,2)
%     plot(Det.Time(Ic),Det.Ang(Ic,2)-90,'.','MarkerSize',ms)
%     
%     h=gca;
%     set(h,'XMinorTick','on')
%     set(h,'XTickLabel',[])
%     set(h,'FontSize',fs);
%     ylabel('Elevation [deg]')
%     grid minor
%     v=axis;
% %     axis([v(1) v(2) 0 90])
%     datetick('x',13)
%     xlabel('Time')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn on brush in figure window, click bad locs while holding down shift
    % key, when done, click space or other kb key
    %
    pause on
    disp('Select Brushing Tool, Brush, Deselect')
    disp('Press Enter to Brush again, or ''q'' to quit')
    pause  % need a keystroke to get out of brush function
    % get key stroke
    cc = get(hf,'CurrentCharacter');
    if strcmp(cc,'q')
        disp('quit')
        break
    end
    
    % get brushed data - 2013b and earlier style
    %     hBrushLine = findall(gca,'tag','Brushing')
    %     brushedData = get(hBrushLine, {'Xdata'});
    %     brushColor = get(hBrushLine, {'Color'});
    
    % get brushed data logical vector - 2016b and later style
    hl = get(gco,'Children');
    if ~isempty(hl)
        x = hl.BrushData;
        brushedIdx = find(hl.BrushData > 0);
    else
        disp('Nothing Brushed, try again')
        continue
    end

    % previous (initial) loop indices
    Ibo = Ib;
    Ico = Ic;
    
    % convert logical vector into indices via find
    if ~isempty(brushedIdx)
        %         brushedIdx = ~isnan(brushedData{1,1});
        Ibi = find(brushedIdx == 1); % brushed = losers
        %aI = find(brushedIdx == 0);  % not brushed = keepers
        Ib = sort([Ibo Ico(Ibi)]);  % combine new brushed and old brush indices
                                    % new brushed values are relative to current plot
        Ic = setdiff(Ia,Ib);    % unbrushed keepers
        disp(['Number of data points kept = ',num2str(length(Ic)),...
            ' out of ',num2str(length(Det.Time)),' locations'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save

% put Indices into structure save
DI.A = Ia;
DI.B = Ib;
DI.C = Ic;

save(sfile,'Det','DI','-mat')



