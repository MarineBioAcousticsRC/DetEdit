function figure50

global dPARAMS dHANDLES p

figure(dHANDLES.spectrafig);clf
dHANDLES.h50 = gca;

if ~isempty(dPARAMS.trueTimes)
    % plot average true click spectrum
    plot(dHANDLES.h50,dPARAMS.ft,dPARAMS.trueSpec,'Linewidth',4)
end

if dPARAMS.ff2 % average false click spec
    % plot average false click spectrum
    hold(dHANDLES.h50, 'on')
    plot(dHANDLES.h50,dPARAMS.ft,dPARAMS.SPEC2,'r','Linewidth',4)
    hold(dHANDLES.h50, 'off')
end
if dPARAMS.ff3  % average id click spec
    hold(dHANDLES.h50, 'on')
       
    for iC = 1:length(dPARAMS.specIDs) % set colors
         dHANDLES.SpecID50{iC} = plot(dHANDLES.h50,dPARAMS.ft,...
             dPARAMS.specID_norm(iC,:),'Linewidth',4);

        set(dHANDLES.SpecID50{iC},'Color',p.colorTab(dPARAMS.specIDs(iC),:),...
            'Visible',dPARAMS.ID_Toggle{iC})
    end
    hold(dHANDLES.h50, 'off')
end


% add figure labels
xlabel(dHANDLES.h50,'Frequency (kHz)');
ylabel(dHANDLES.h50,'Normalized Received Level')
grid(dHANDLES.h50,'on')
xlim(dHANDLES.h50, 'manual');
ylim(dHANDLES.h50,[0 1]);
xlim(dHANDLES.h50,[p.fLow,p.fHi])

if ~isempty(dPARAMS.yell) && ~isempty(dPARAMS.csnJ)
    
    hold(dHANDLES.h50,'on') % add click to spec plot in BLACK
    cspJy = mean(dPARAMS.cspJ(dPARAMS.yell,:),1);
    tSPEC = norm_spec_simple(cspJy,dPARAMS.fimint,dPARAMS.fimaxt);
    plot(dHANDLES.h50,dPARAMS.ft,tSPEC,'k','Linewidth',4);
    hold(dHANDLES.h50,'off')
end