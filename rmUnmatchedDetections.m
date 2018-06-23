function detLabels = rmUnmatchedDetections(MTT, detLabels)

jDet = []; ia = []; ic = [];
if (~isempty(detLabels))
    [jDet,~,~] = intersect(MTT,detLabels);
    rDet = length(detLabels) - length(jDet);
    disp([' Removed ',num2str(rDet),...
        ' detections that do not match detection times']);
    if size(jDet,1)<size(jDet,2)
        jDet = jDet';
    end
    detLabels = jDet;
end