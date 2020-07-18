function ID_toggle(hObject,eventdata,iColor)

global dPARAMS p dHANDLES
dPARAMS.ID_Toggle(iColor) = get(dHANDLES.h10handles.colorbox{iColor},'Value');
if dPARAMS.ID_Toggle
    stateString = 'ON';
else 
    stateString = 'OFF';
end

fprintf('Toggled %s %s\n',get(dHANDLES.h10handles.spLabel{iColor},'String'),...
    stateString)


figure201
figure51
figure53