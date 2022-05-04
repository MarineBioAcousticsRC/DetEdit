function figure50

global dPARAMS dHANDLES p

figure(dHANDLES.spectrafig);clf
dHANDLES.h50 = gca;

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