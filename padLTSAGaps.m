function [ pt, pwr ] = padLTSAGaps(hdr,L,pwr)
%
% Takes LTSA raw file start times and returns a vector of start times for
% the individual spectral averages.  Accounts for duty cycle
%
% pt = vector containing start times of spectral averages
%
% hdr = Struct containing LTSA metadata
% L = vector of raw file indexes
%

pt = [];
mnum2secs = 24*60*60;
Y2K = datenum([ 2000 0 0 0 0 0 ]);
date_fmt = 'mm/dd/yy HH:MM:SS';

if hdr.ltsa.dnumStart(L)< Y2K
    rfStarts = hdr.ltsa.dnumStart(L) + Y2K;
else
    rfStarts = hdr.ltsa.dnumStart(L);
end
diffRF_s = round(diff(rfStarts)*mnum2secs);
nbin = sum(hdr.ltsa.nave(L));
t1 = rfStarts(1);
dt_dnum = datenum([ 0 0 0 0 0 hdr.ltsa.tave ]);
% Assume max raw file size in set is the standard. Anything less should be padded.
rfdt = hdr.ltsa.dur(L);%rfdt = max((hdr.ltsa.dnumEnd - hdr.ltsa.dnumStart)*(24*60*60));

gapsidx = find(diffRF_s~=rfdt(1:end-1));
rfdt1_dnum = datenum([ 0 0 0 0 0 rfdt(1)-5 ]);
if size(unique(diffRF_s), 2) > 1
    fprintf('Duty cycle or gap encountered!\n')
    fprintf('\tBout starting at %s\n', datestr(rfStarts(1),date_fmt));
    pt = t1:dt_dnum:t1+rfdt1_dnum; % make first raw file's part of pt;
    rf = 2;
    while rf <= length(L)
        if ismember(rf,gapsidx+1) % if there's a duty cycle or gap, pad pwr for later
            tdt = round((rfStarts(rf) - rfStarts(rf-1))*mnum2secs-rfdt(rf));
            numaves2pad = ceil(tdt/hdr.ltsa.tave);
            pwrpad = ones(size(pwr,1),numaves2pad).*-128;
            gaveidx = length(pt); % last average before gap
            try
                pwr = [ pwr(:,1:gaveidx), pwrpad, pwr(:,gaveidx+1:end)];
            catch
                warning('unknown error padding LTSA time gap, ouput may be incorrect')
                rf = rf+1;
                continue
            end
            padavetimes = pt(end)+dt_dnum:dt_dnum:rfStarts(rf)-dt_dnum;
            pt = [ pt, padavetimes ];
        end
        rfdt_dnum = datenum([ 0 0 0 0 0 rfdt(rf)-5 ]);
        tpt = rfStarts(rf):dt_dnum:rfStarts(rf)+rfdt_dnum;
        pt = [ pt, tpt ];
        rf = rf+1;
    end
    if size(pt,2)>size(pwr,2)
        % if size still doesn't match add zeros to pwr.        
        pwr = [pwr,zeros(size(pwr,1),size(pt,2)-size(pwr,2))];

    end
else
    pt = [t1:dt_dnum:t1 + (nbin-1)*dt_dnum];
end

end