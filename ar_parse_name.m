function info = ar_parse_name(filename)

% ar_parse_name
%
% Split a TPWS / ID / LTSA file name into deployment and disk, so the three
% kinds of files can be matched even when the rest of their names differ.
%
%   NFC_A_02_disk01a_combinedBWD_TPWS1.mat -> dep 'nfca02', disk 1, letter 'a'
%   NFC_A_02_disk01a_combinedBWD_ID1.mat   -> dep 'nfca02', disk 1, letter 'a'
%   NFC_A_02_disk01_5s_100Hz.ltsa          -> dep 'nfca02', disk 1, letter ''
%   NFC_A_05_disk04a1_combinedBWD_TPWS1    -> dep 'nfca05', disk 4, letter 'a1'
%     (so disk04a and disk04a1 are different disks)
%
% The deployment is everything before "disk", lowercased with separators
% (_ - . spaces) removed, so NFC_A_02 and NFC-A-02 match.
%
% Output struct:
%   ok      - false if the name has no diskNN tag
%   depRaw  - deployment as written (NFC_A_02)
%   dep     - normalized deployment (nfca02)
%   diskNum - disk number (1)
%   letter  - lowercase sub-disk tag ('a' or 'a1'), '' if none
%   key     - dep + disk + letter, e.g. 'nfca02|disk01a'

[~,stem] = fileparts(filename);
info = struct('ok',false,'depRaw','','dep','','diskNum',NaN,'letter','','key','');
tok = regexp(stem,'^(.*?)[^a-zA-Z0-9]*[dD][iI][sS][kK](\d+)([a-zA-Z]\d*(?![a-zA-Z0-9]))?',...
    'tokens','once');
if isempty(tok) || isempty(tok{1})
    return
end
info.ok = true;
info.depRaw = tok{1};
info.dep = lower(regexprep(tok{1},'[^a-zA-Z0-9]',''));
info.diskNum = str2double(tok{2});
if numel(tok) >= 3
    info.letter = lower(tok{3});
end
info.key = sprintf('%s|disk%02d%s',info.dep,info.diskNum,info.letter);
