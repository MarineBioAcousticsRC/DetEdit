function [tf,tffreq,tfuppc] = getTransfunc(filePrefix, tfName,p)

if exist('tfName','var')
    disp('Load Transfer Function File');
    stndeploy = strsplit(filePrefix,'_'); % get only station and deployment
    tffn = findTfFile(tfName,stndeploy); % get corresponding tf file
else
    error('No path specified for transfer function files. Add tfName')
end

fid = fopen(tffn);
[A,count] = fscanf(fid,'%f %f',[2,inf]);
tffreq = A(1,:);
tfuppc = A(2,:);
fclose(fid);

tf = interp1(tffreq,tfuppc,p.tfSelect,'linear','extrap');
disp(['TF @',num2str(p.tfSelect),' Hz =',num2str(tf)]);
