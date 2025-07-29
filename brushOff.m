function brushOff(varargin)
% when mouse is released, turn off 
% brush so we can listen for keyboard command

% also store state, so that we can turn 
% brush back on after command is recieved.
figH = varargin{3};
hManager = uigetmodemanager(figH);
set(hManager.CurrentMode,'KeyPressFcn',@keyAction)
%%% solution from:
%%% https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
% Overrides default behavior to allow listening
% for key press with brush enabled.





% global dHANDLES
% figName = varargin{4};
% 
% brushInfo = brush(figH);
% if strcmp(figName,'LTSA')
%     % record state of brush for later.
%     dHANDLES.LTSA.brushState = brushInfo.Enable;
%     dHANDLES.LTSA.brushColor = brushInfo.Color;
%     disp('brush state stored.')
%     dHANDLES.LTSA.brushState
% end
% figure(figH)
% brush off


