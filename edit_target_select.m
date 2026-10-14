function edit_target_select(hObject,~)

% edit_target_select
%
% Callback for the "Labels edit:" menu in the ID legend. Chooses what hand
% labeling (brush colors, number keys, ALT + digits, 'l', 'e') changes:
%   Model + GT - model output (col 2) and ground truth (col 4)
%   Model only - model output only
%   GT only    - ground truth only
% See set_ID_labels.m.

global dPARAMS dHANDLES

targets = {'both','model','gt'};
names = {'model output + ground truth','model output only','ground truth only'};
iT = get(hObject,'Value');
dPARAMS.editTarget = targets{iT};
fprintf('Hand labels now edit: %s\n',names{iT});
boutMotion % redraw so the LTSA title shows the current choice

% give keyboard focus back to the plots so key commands keep working
set(hObject,'Enable','off');
drawnow;
set(hObject,'Enable','on');
if isfield(dHANDLES,'LTSAfig') && ishghandle(dHANDLES.LTSAfig)
    figure(dHANDLES.LTSAfig)
end
