function figure52

global dPARAMS dHANDLES p

figure(dHANDLES.wavefig)
clf;
dHANDLES.h52 = gca;

if p.specploton
    hold(dHANDLES.h52, 'on')
    if ~isempty(dPARAMS.trueTimes) % average true waveform in blue
        dHANDLES.WAV52 = plot(dHANDLES.h52, dPARAMS.wtrue,'Visible',dPARAMS.NoLabel_Toggle);
    end
    
    dHANDLES.WAVFD52 = [];
    if dPARAMS.ff2   % average false waveform in red
        dHANDLES.WAVFD52 = plot(dHANDLES.h52,dPARAMS.wavFD + 0.5 ,'r',...
            'Visible',dPARAMS.FD_Toggle);
    end
    
    dHANDLES.WAVID52 = cell(size(p.colorTab,1),1);
    if dPARAMS.ff3 % average waveform ID'd in associated color
        for iC = 1:size(dPARAMS.wavID,1) % set colors
            iColor = dPARAMS.specIDs(iC);
            dHANDLES.WAVID52{iColor} = plot(dHANDLES.h52,(dPARAMS.wavID(iC,:)' +...,
                -ones(size(dPARAMS.wavID(iC,:),1),1)'.*rand(1,min(size(dPARAMS.wavID(iC,:),1)))));
            set(dHANDLES.WAVID52{iColor},'Color',p.colorTab(dPARAMS.specIDs(iC),:),...
                'Visible',dPARAMS.ID_Toggle{iColor})
        end
    end
    hold(dHANDLES.h52, 'off')
end

xlabel(dHANDLES.h52,sprintf('Samples (%0.1fms @ %0.0fkHz)',...
    (size(dPARAMS.wtrue,2)/p.sampleRate),p.sampleRate));
ylabel(dHANDLES.h52,' Normalized Amplitude');

% if you have items brushed in yellow, highlight those on each plot
if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h52,'on') % add click to waveform plot in BLACK
    plot(dHANDLES.h52,norm_wav(mean(dPARAMS.csnJ(dPARAMS.yell,:),1))' + 1.5,'k');
    hold(dHANDLES.h52,'off')
    
end