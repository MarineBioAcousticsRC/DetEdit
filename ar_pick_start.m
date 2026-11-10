function idx = ar_pick_start(files, filePrefix)

% ar_pick_start
%
% Accept/reject review: which disk detEdit opens first, from filePrefix.
%   'NFC_A_02_disk01a' (or any name with that deployment + disk) -> that disk
%   'NFC_A_02_disk01'  -> the first sub-disk of disk01 (disk01, disk01a, ...)
%   'NFC_A_02'         -> the first disk of that deployment
%   ''                 -> the first disk overall
% Deployments are compared like ar_parse_name does (case and _ - . ignored).
%
% Inputs:
%   files - struct array from ar_build_index
%   filePrefix - p.filePrefix
% Output:
%   idx - index into files

if isempty(files)
    error('detEdit:acceptReject',['No disks to review: no accepted/rejected ',...
        'ID file matched a TPWS file.']);
end
idx = 1;
if isempty(filePrefix)
    return
end
info = ar_parse_name(filePrefix);
if info.ok
    hit = find(strcmp({files.dep},info.dep) & [files.diskNum] == info.diskNum);
    if ~isempty(info.letter)
        hit = hit(strcmp({files(hit).letter},info.letter));
    end
else
    dep = lower(regexprep(filePrefix,'[^a-zA-Z0-9]',''));
    hit = find(strcmp({files.dep},dep));
end
if isempty(hit)
    error('detEdit:acceptReject',['filePrefix ''%s'' does not match any disk with ',...
        'accepted/rejected labels. First disk is %s.'],filePrefix,files(1).label);
end
idx = hit(1);
