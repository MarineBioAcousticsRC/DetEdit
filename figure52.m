function figure52

global dPARAMS dHANDLES p

figure(dHANDLES.wavefig)
clf;
dHANDLES.h52 = gca;

if p.specploton
    if ~isempty(dPARAMS.trueTimes)
        % average true click waveform
        plot(dHANDLES.h52, dPARAMS.wtrue);
    end
    if dPARAMS.ff2   % average false click spec
        % plot average false click waveform
        hold(dHANDLES.h52, 'on')
        plot(dHANDLES.h52,dPARAMS.wavFD + 0.5 ,'r');
        hold(dHANDLES.h52, 'off')
    end
    if dPARAMS.ff3
        % plot average ID'd click waveform(s)
        hold(dHANDLES.h52, 'on')
        hID2 = plot(dHANDLES.h52,(wavID + repmat(-1*rand(size(hID)),1,size(wavID,2)))');
        
        for iC = 1:length(hID) % set colors
            set(hID2(iC),'Color',p.colorTab(dPARAMS.specIDs(iC),:))
        end
        hold(dHANDLES.h52, 'off')
    end
end

xlabel(dHANDLES.h52,sprintf('Samples (%0.1fms @ %0.0fkHz)',...
    (size(dPARAMS.wtrue,2)/p.sampleRate),p.sampleRate));
ylabel(dHANDLES.h52,' Normalized Amplitude');

% if you have items brushed in yellow, highlight those on each plot
if p.specploton && ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h52,'on') % add click to waveform plot in BLACK
    plot(dHANDLES.h52,dPARAMS.norm_wav(mean(dPARAMS.csnJ(dPARAMS.yell,:),1))' + 1.5,'k');
    hold(dHANDLES.h52,'off')
    
end