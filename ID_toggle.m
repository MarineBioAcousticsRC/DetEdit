function ID_toggle(hObject,eventdata,iColor)

global dPARAMS dHANDLES

if iColor == 99
    % toggle falses
    stateString = get(dHANDLES.h10handles.colorboxFD,'Value');
    if stateString
        dPARAMS.FD_Toggle = 'on';
    else
        dPARAMS.FD_Toggle = 'off';
    end
    set(dHANDLES.RLFD201,'Visible',dPARAMS.FD_Toggle)
    set(dHANDLES.ICIFD201,'Visible',dPARAMS.FD_Toggle)
    set(dHANDLES.RLFD51,'Visible',dPARAMS.FD_Toggle)
    set(dHANDLES.RMSFD53,'Visible',dPARAMS.FD_Toggle)
    fprintf('Toggled False %s\n',dPARAMS.FD_Toggle)
elseif iColor == 0
    % toggle unlabeled
    stateString = get(dHANDLES.h10handles.colorboxUnlabeled,'Value');
    if stateString
        dPARAMS.NoLabel_Toggle = 'on';
    else
        dPARAMS.NoLabel_Toggle = 'off';
    end
    set(dHANDLES.RL201,'Visible',dPARAMS.NoLabel_Toggle)
    set(dHANDLES.ICI201,'Visible',dPARAMS.NoLabel_Toggle)
    set(dHANDLES.RL51,'Visible',dPARAMS.NoLabel_Toggle)
    set(dHANDLES.RMS53,'Visible',dPARAMS.NoLabel_Toggle)
    
    fprintf('Toggled Unlabeled %s\n',dPARAMS.NoLabel_Toggle)
else
    % toggle ID
    stateString = get(dHANDLES.h10handles.colorbox{iColor},'Value');
    if iColor<=size(dHANDLES.RLID201,2)
        
        if stateString
            dPARAMS.ID_Toggle{iColor}  = 'on';
        else
            dPARAMS.ID_Toggle{iColor}  = 'off';
        end
         set(dHANDLES.RLID201{iColor},'Visible',dPARAMS.ID_Toggle{iColor})
        set(dHANDLES.ICIID201{iColor},'Visible',dPARAMS.ID_Toggle{iColor})
        set(dHANDLES.RLID51{iColor},'Visible',dPARAMS.ID_Toggle{iColor})
        set(dHANDLES.RMSID53{iColor},'Visible',dPARAMS.ID_Toggle{iColor})
        fprintf('Toggled %s %s\n',get(dHANDLES.h10handles.spLabel{iColor},'String'),...
            dPARAMS.ID_Toggle{iColor})
    end
end