function [nb,eb,sb,bd] = ar_calculate_bouts(boutTimes,p)

% ar_calculate_bouts
%
% Accept/reject review: sessions from the labeled click times
% (ar_bout_times). Used by both mkLTSAsessions and detEdit so the session
% lists match.
%
% Differences from calculate_bouts:
%   - minBout = 0 keeps every session, even one made of a single labeled
%     click (calculate_bouts drops zero-length bouts), so no reviewed click
%     is left without a session
%   - sessions shorter than p.arMinSession minutes are widened to that
%     length (centered), so you see the clicks around them and the LTSA
%     panel has a time span to draw
%
% Outputs as in calculate_bouts.

q = p;
if isempty(q.minBout) || q.minBout == 0
    q.minBout = [];
end
[nb,eb,sb,bd] = calculate_bouts(boutTimes,q);

minLen = 10/(24*60); % default 10 minutes
if isfield(p,'arMinSession') && ~isempty(p.arMinSession)
    minLen = p.arMinSession/(24*60);
end
short = bd < minLen;
mid = (sb(short) + eb(short))/2;
sb(short) = mid - minLen/2;
eb(short) = mid + minLen/2;
bd = eb - sb;
