function SPEC = norm_spec(spectraSet, flowt,fimint,fimaxt)
% normalize spectrum given cutoff frequencies
mnSPEC = min(spectraSet(:,fimint+5:flowt),[],2);
SPEC = spectraSet - mnSPEC;  % make low freq part = 0
mxSPEC = max(SPEC(:,flowt:fimaxt),[],2);
SPEC = SPEC ./ mxSPEC;  % make high freq part = 1


