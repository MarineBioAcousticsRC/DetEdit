function [y] = groupoverlaps(x)
% Nx2 matrix of endpoints x1, x2 of intervals

% -----plot it first
% figure
% nline  = repmat((1:size(x,1))',1,2);
% plot(x',nline','o-')

% ----- find union of intervals
x = sort(x,2); % order start and end
nrow = size(x,1);
[x, ind] = sort(x(:));
n = [(1:nrow) (1:nrow)]';
n = n(ind);
c = [ones(1,nrow) -ones(1,nrow)]';
c = c(ind);
csc = cumsum(c);          % =0 at upper end of new interval(s)
irit = find(csc==0);
ilef = [1; irit+1];
ilef(end) = [];           % no new interval starting at the very end
% y matrix is start and end points of the new intervals, y1,y2
%
% ny matrix is the corresponding indices of the start and end points
% in terms of what row of x they occurred in.
y = [x(ilef) x(irit)];
% ny = [n(ilef) n(irit)];
% 
% figure(2)
% nline  = repmat((1:size(y,1))',1,2);
% plot(y',nline','o-')



  % problem with 20 - GofMX_MC01_disk01-16.mat and 17 - GofMX_GC06_disk01-15.mat