function keyRelease(varargin)

% keyRelease
%
% Finishes an ALT + digits label code. While ALT is held, keyAction collects
% the digits typed (dPARAMS.altDigits). When ALT is released, the collected
% code (e.g. '12') is applied to the brushed points exactly like pressing a
% single number key does for labels 1-9.

global dPARAMS

evtKey = '';
if nargin >= 2 && ~isempty(varargin{2})
    try
        evtKey = varargin{2}.Key;
    catch
    end
end
if ~strcmp(evtKey,'alt')
    return
end
if ~isfield(dPARAMS,'altDigits') || isempty(dPARAMS.altDigits)
    return % ALT pressed and released without digits: nothing to do
end

code = str2double(dPARAMS.altDigits);
dPARAMS.altDigits = '';
fprintf('Applying label %d to brushed points\n',code);

% hand the code to keyAction as if it had been typed as one key
dPARAMS.pendingCC = num2str(code);
keyAction(varargin{1},[])
