function brushOff(varargin)
% brushOff
%
% Brush ActionPostCallback: when the mouse is released after brushing,
% make sure DetEdit's key commands (and ALT + digits) still reach keyAction
% / keyRelease while brush mode stays on.
%
% Brush mode normally swallows key presses. Two approaches are tried, so
% this works in older MATLAB (e.g. R2016b) and newer releases (R2024a/b):
%   1) give the active mode our key handlers (what worked in R2016b)
%   2) switch off the mode's key listeners and use the figure's own
%      handlers (https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan)
% Each step is wrapped so a change in MATLAB's internals can't stop DetEdit.

figH = varargin{3};

try
    hManager = uigetmodemanager(figH);
catch
    return % mode manager not available: nothing more we can do
end

% 1) older approach: key handlers on the active mode
try
    set(hManager.CurrentMode,'KeyPressFcn',@keyAction)
catch
end
try
    set(hManager.CurrentMode,'KeyReleaseFcn',@keyRelease)
catch
end

% 2) newer approach: disable the mode's window key listeners, then the
% figure's own KeyPressFcn / KeyReleaseFcn are used
try
    [hManager.WindowListenerHandles.Enabled] = deal(false);
catch
    try
        set(hManager.WindowListenerHandles,'Enable','off'); % very old MATLAB
    catch
    end
end
try
    set(figH,'WindowKeyPressFcn',[]);
catch
end
try
    set(figH,'KeyPressFcn',@keyAction,'KeyReleaseFcn',@keyRelease);
catch
end
