function stack_ID_rings(ax)

% stack_ID_rings
%
% Reorder the ID circles drawn by plot_misclass_rings in one axes so that
% mistakes are never hidden under correct points:
%   top    - wrong (ground truth differs from the prediction)
%   middle - no ground truth
%   bottom - correct (ground truth matches the prediction), and everything
%            else in the axes keeps its existing order beneath those
%
% Call after all ID points for the axes have been plotted.

ch = get(ax,'Children');
if isempty(ch)
    return
end
gtMatch = cell(size(ch));
for iC = 1:length(ch)
    if isappdata(ch(iC),'gtMatch')
        gtMatch{iC} = getappdata(ch(iC),'gtMatch');
    else
        gtMatch{iC} = '';
    end
end
isWrong = strcmp(gtMatch,'wrong');
isNone = strcmp(gtMatch,'none');
isCorrect = strcmp(gtMatch,'correct');
isOther = ~isWrong & ~isNone & ~isCorrect;

% first child is drawn on top; correct ones sit just above everything else
set(ax,'Children',[ch(isWrong); ch(isNone); ch(isCorrect); ch(isOther)]);
