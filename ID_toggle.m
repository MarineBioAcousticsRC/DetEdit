function ID_toggle(hObject,eventdata,iColor)

global dPARAMS dHANDLES

if iColor == 99
    % toggle falses
    dPARAMS.FD_Toggle = get(dHANDLES.h10handles.colorboxFD,'Value');
    if dPARAMS.FD_Toggle
        set(dHANDLES.RLFD201,'Visible',1)
        set(dHANDLES.ICIFD201,'Visible',1)
        set(dHANDLES.RLFD51,'Visible',1)
        set(dHANDLES.RMSFD53,'Visible',1)
        stateString = 'ON';
    else
        set(dHANDLES.RLFD201,'Visible',0)
        set(dHANDLES.ICIFD201,'Visible',0)
        set(dHANDLES.RLFD51,'Visible',0)
        set(dHANDLES.RMSFD53,'Visible',0)
        stateString = 'OFF';
    end
    fprintf('Toggled False %s\n',stateString)
elseif iColor == 0
    % toggle unlabeled
    dPARAMS.NoLabel_Toggle = get(dHANDLES.h10handles.colorboxUnlabeled,'Value');
    if dPARAMS.NoLabel_Toggle
        set(dHANDLES.RL201,'Visible',1)
        set(dHANDLES.ICI201,'Visible',1)
        set(dHANDLES.RL51,'Visible',1)
        set(dHANDLES.RMS53,'Visible',1)
        stateString = 'ON';
    else
        set(dHANDLES.RL201,'Visible',0)
        set(dHANDLES.ICI201,'Visible',0)
        set(dHANDLES.RL51,'Visible',0)
        set(dHANDLES.RMS53,'Visible',0)
        stateString = 'OFF';
    end
    fprintf('Toggled Unlabeled %s\n',stateString)
else
    % toggle ID
    dPARAMS.ID_Toggle(iColor) = get(dHANDLES.h10handles.colorbox{iColor},'Value');
    
    if iColor<=size(dHANDLES.RLID201,2)
        
        if dPARAMS.ID_Toggle
            set(dHANDLES.RLID201{iColor},'Visible',1)
            set(dHANDLES.ICIID201{iColor},'Visible',1)
            set(dHANDLES.RLID51{iColor},'Visible',1)
            set(dHANDLES.RMSID53{iColor},'Visible',1)
            stateString = 'ON';
        else
            set(dHANDLES.RLID201{iColor},'Visible',0)
            set(dHANDLES.ICIID201{iColor},'Visible',0)
            set(dHANDLES.RLID51{iColor},'Visible',0)
            set(dHANDLES.RMSID53{iColor},'Visible',0)
            stateString = 'OFF';
        end
        fprintf('Toggled %s %s\n',get(dHANDLES.h10handles.spLabel{iColor},'String'),...
            stateString)
    end
end