function figure52

global dPARAMS dHANDLES p

% With p.plotCorrect_incorrect on, this window has two panels:
%   top    - average waveforms by model output (predicted label, zID col 2),
%            same as the single-panel version
%   bottom - average waveforms by ground truth label (zID col 4), each on
%            its own row so they don't overlap
% The window is made twice as tall the first time, so each panel keeps
% about the shape of the original single plot.
showGTpanel = isfield(p,'plotCorrect_incorrect') && p.plotCorrect_incorrect;
waveSpacing = 1.2; % vertical offset between stacked ground truth waveforms
waveLW = 1.5; % waveform line width (MATLAB default is 0.5)

figure(dHANDLES.wavefig)
clf;
if showGTpanel
    if ~isfield(dHANDLES,'waveTall') || ~dHANDLES.waveTall
        make_window_taller(dHANDLES.wavefig)
        dHANDLES.waveTall = true;
    end
    dHANDLES.h52 = axes('Parent',dHANDLES.wavefig,'Position',[.13 .57 .775 .37]);
    dHANDLES.h52gt = axes('Parent',dHANDLES.wavefig,'Position',[.13 .08 .775 .37]);
else
    dHANDLES.h52 = gca;
    dHANDLES.h52gt = [];
end

if p.specploton
    hold(dHANDLES.h52, 'on')
    if ~isempty(dPARAMS.trueTimes) % average true waveform in blue
        dHANDLES.WAV52 = plot(dHANDLES.h52, dPARAMS.wtrue,'Visible',dPARAMS.NoLabel_Toggle,...
            'LineWidth',waveLW);
    end

    dHANDLES.WAVFD52 = [];
    if dPARAMS.ff2   % average false waveform in red
        dHANDLES.WAVFD52 = plot(dHANDLES.h52,dPARAMS.wavFD + 0.5 ,'r',...
            'Visible',dPARAMS.FD_Toggle,'LineWidth',waveLW);
    end

    dHANDLES.WAVID52 = cell(size(p.colorTab,1),1);
    if dPARAMS.ff3 % average waveform ID'd in associated color
        for iC = 1:size(dPARAMS.wavID,1) % set colors
            iColor = dPARAMS.specIDs(iC);
            dHANDLES.WAVID52{iColor} = plot(dHANDLES.h52,(dPARAMS.wavID(iC,:)' +...,
                -ones(size(dPARAMS.wavID(iC,:),1),1)'.*rand(1,min(size(dPARAMS.wavID(iC,:),1)))));
            set(dHANDLES.WAVID52{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC),:),...
                'Visible',dPARAMS.ID_Toggle{iColor},'LineWidth',waveLW)
        end
    end
    hold(dHANDLES.h52, 'off')

    % bottom panel: average waveforms by ground truth label, one row each
    nGT = 0;
    if showGTpanel
        hold(dHANDLES.h52gt,'on')
        if dPARAMS.ff3 && ~isempty(dPARAMS.gtCode)
            % same clicks as the top panel's ID lines: above the confidence limit
            gtSet = dPARAMS.gtCode;
            gtSet(~dPARAMS.labelConfIdx) = NaN;
            gtIDs = unique(gtSet(~isnan(gtSet) & gtSet <= size(p.colorTab,1)));
            nGT = length(gtIDs);
            for iG = 1:nGT
                thisGT = gtIDs(iG);
                rows = dPARAMS.K3(gtSet == thisGT);
                gtWav = norm_wav(mean(dPARAMS.csnJ(rows,:),1)) + waveSpacing*(nGT - iG);
                hGT = plot(dHANDLES.h52gt,gtWav,'Color',p.colorTab(thisGT,:),'LineWidth',waveLW,...
                    'Visible',dPARAMS.ID_Toggle{thisGT});
                % legend toggle for this label hides its line in both panels
                dHANDLES.WAVID52{thisGT} = [dHANDLES.WAVID52{thisGT}; hGT];
            end
        end
        hold(dHANDLES.h52gt,'off')
        if nGT == 0
            title(dHANDLES.h52gt,'Ground truth (none in this session)')
        else
            title(dHANDLES.h52gt,'Ground truth')
            ylim(dHANDLES.h52gt,[-.05, waveSpacing*nGT + 1.05]) % room for highlight row
            set(dHANDLES.h52gt,'YTick',[]) % stacked rows: amplitude ticks don't apply
        end
        title(dHANDLES.h52,'Model output')
    end
end

for ax = [dHANDLES.h52, dHANDLES.h52gt]
    xlabel(ax,sprintf('Samples (%0.1fms @ %0.0fkHz)',...
        (size(dPARAMS.wtrue,2)/p.sampleRate),p.sampleRate));
    ylabel(ax,' Normalized Amplitude');
end

% if you have items brushed in yellow, highlight those on each plot
if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    yellWav = norm_wav(mean(dPARAMS.csnJ(dPARAMS.yell,:),1))';
    hold(dHANDLES.h52,'on') % add click to waveform plot in BLACK
    plot(dHANDLES.h52,yellWav + 1.5,'k','LineWidth',waveLW);
    hold(dHANDLES.h52,'off')
    if ~isempty(dHANDLES.h52gt)
        hold(dHANDLES.h52gt,'on') % and on its own top row in the ground truth panel
        plot(dHANDLES.h52gt,yellWav + waveSpacing*nGT,'k','LineWidth',waveLW);
        hold(dHANDLES.h52gt,'off')
    end
end
