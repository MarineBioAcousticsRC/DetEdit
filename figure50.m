function figure50

global dPARAMS dHANDLES p

% With p.plotCorrect_incorrect on, this window has two panels:
%   top    - average spectra by model output (predicted label, zID col 2),
%            same as the single-panel version
%   bottom - average spectra by ground truth label (zID col 4)
% The window is made twice as tall the first time, so each panel keeps
% about the shape of the original single plot.
showGTpanel = isfield(p,'plotCorrect_incorrect') && p.plotCorrect_incorrect;

figure(dHANDLES.spectrafig);clf
if showGTpanel
    if ~isfield(dHANDLES,'spectraTall') || ~dHANDLES.spectraTall
        make_window_taller(dHANDLES.spectrafig)
        dHANDLES.spectraTall = true;
    end
    dHANDLES.h50 = axes('Parent',dHANDLES.spectrafig,'Position',[.13 .57 .775 .37]);
    dHANDLES.h50gt = axes('Parent',dHANDLES.spectrafig,'Position',[.13 .08 .775 .37]);
else
    dHANDLES.h50 = gca;
    dHANDLES.h50gt = [];
end

if p.specploton
    hold(dHANDLES.h50, 'on')
    if ~isempty(dPARAMS.trueTimes)% average true spec in blue
        dHANDLES.SPE50 = plot(dHANDLES.h50,dPARAMS.ft,dPARAMS.trueSpec,'Linewidth',4,...
            'Visible',dPARAMS.NoLabel_Toggle);
    end

    dHANDLES.SPEFD50 = [];
    if dPARAMS.ff2 % average false spec in red
        dHANDLES.SPEFD50 = plot(dHANDLES.h50,dPARAMS.ft,dPARAMS.SPEC2,'r',...
            'Linewidth',4,'Visible',dPARAMS.FD_Toggle);
    end

    dHANDLES.SPEID50 = cell(size(p.colorTab,1),1);
    if dPARAMS.ff3  % average id click spec
        for iC = 1:length(dPARAMS.specIDs) % set colors
            iColor = dPARAMS.specIDs(iC);
            dHANDLES.SPEID50{iColor} = plot(dHANDLES.h50,dPARAMS.ft,...
                dPARAMS.specID_norm(iC,:),'Linewidth',4);
            set(dHANDLES.SPEID50{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC),:),...
                'Visible',dPARAMS.ID_Toggle{iColor})
        end
    end
    hold(dHANDLES.h50, 'off')

    % bottom panel: average spectra by ground truth label
    if showGTpanel
        nGT = 0;
        hold(dHANDLES.h50gt,'on')
        if dPARAMS.ff3 && ~isempty(dPARAMS.gtCode)
            % same clicks as the top panel's ID lines: above the confidence limit
            gtSet = dPARAMS.gtCode;
            gtSet(~dPARAMS.labelConfIdx) = NaN;
            gtIDs = unique(gtSet(~isnan(gtSet) & gtSet <= size(p.colorTab,1)));
            for iG = 1:length(gtIDs)
                thisGT = gtIDs(iG);
                rows = dPARAMS.K3(gtSet == thisGT);
                gtSpec = norm_spec_simple(mean(dPARAMS.cspJ(rows,:),1),...
                    dPARAMS.fimint,dPARAMS.fimaxt);
                hGT = plot(dHANDLES.h50gt,dPARAMS.ft,gtSpec,'Linewidth',4,...
                    'Color',p.colorTab(thisGT,:),'Visible',dPARAMS.ID_Toggle{thisGT});
                % legend toggle for this label hides its line in both panels
                dHANDLES.SPEID50{thisGT} = [dHANDLES.SPEID50{thisGT}; hGT];
                nGT = nGT + 1;
            end
        end
        hold(dHANDLES.h50gt,'off')
        if nGT == 0
            title(dHANDLES.h50gt,'Ground truth (none in this session)')
        else
            title(dHANDLES.h50gt,'Ground truth')
        end
        title(dHANDLES.h50,'Model output')
    end
end

% add figure labels
for ax = [dHANDLES.h50, dHANDLES.h50gt]
    xlabel(ax,'Frequency (kHz)');
    ylabel(ax,'Normalized Received Level')
    grid(ax,'on')
    xlim(ax, 'manual');
    ylim(ax,[0 1]);
    xlim(ax,[p.fLow,p.fHi])
end

if ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    % add highlighted clicks' spectrum in BLACK (both panels)
    cspJy = mean(dPARAMS.cspJ(dPARAMS.yell,:),1);
    tSPEC = norm_spec_simple(cspJy,dPARAMS.fimint,dPARAMS.fimaxt);
    for ax = [dHANDLES.h50, dHANDLES.h50gt]
        hold(ax,'on')
        plot(ax,dPARAMS.ft,tSPEC,'k','Linewidth',4);
        hold(ax,'off')
    end
end
end
