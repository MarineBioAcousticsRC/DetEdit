function make_window_taller(fig)

% make_window_taller
%
% Used by figure50/figure52 when the ground truth panel is added.
% double the window's height (keeping its top edge), within the screen
oldUnits = get(fig,'Units');
set(fig,'Units','pixels');
pos = get(fig,'Position');
screen = get(0,'ScreenSize');
newH = 2*pos(4);
if screen(4) > 200 % only cap to the screen if its size is known
    newH = min(newH, screen(4) - 100);
end
newH = max(newH, pos(4)); % never shrink
newY = max(pos(2) + pos(4) - newH, 40);
set(fig,'Position',[pos(1), newY, pos(3), newH]);
set(fig,'Units',oldUnits);
end
