function muSPEC = norm_spec(spectraSet, flowt,fimint,fimaxt)
% normalize spectrum given cutoff frequencies
% normalizes each spectrum independently, then averages,
mnSPEC = min(spectraSet(:,fimint:flowt),[],2);
SPEC = spectraSet - repmat(mnSPEC, 1,size(spectraSet,2));  % make low freq part = 0
mxSPEC = max(SPEC(:,flowt:fimaxt),[],2);
SPEC2 = SPEC ./ repmat(mxSPEC,1,size(spectraSet,2));  % make high freq part = 1
goodRows = ~isnan(mean(SPEC2,2));
muSPEC1 = nanmean(SPEC2(goodRows,:),1);

% then rescale so that all peaks hit 1

muSPEC = (muSPEC1-min(muSPEC1))./max(muSPEC1-min(muSPEC1));


