function muSPEC = norm_spec(spectraSet, flowt,fimint,fimaxt)
% normalize spectrum given cutoff frequencies
mnSPEC = min(spectraSet,[],2);
SPEC = spectraSet - repmat(mnSPEC, 1,size(spectraSet,2));  % make low freq part = 0
mxSPEC = max(SPEC,[],2);
SPEC2 = SPEC ./ repmat(mxSPEC,1,size(spectraSet,2));  % make high freq part = 1
goodRows = ~isnan(mean(SPEC2,2));
muSPEC = nanmean(SPEC2(goodRows,:),1);


