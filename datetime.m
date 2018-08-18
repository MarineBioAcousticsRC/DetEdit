function d = datetime(t);

% DATETIME  returns correct DateTime from Matlab's DATENUM
%
% Takes care on accuraccy-problems of DATEVEC, 
%    which returns seconds == 59.99999999999272
%
% Date = DATETIME( Day )
%
% Day   Matlab's DATENUM
% Date  [ YY MM DD hh mm ss ]   6 Columns
%
% DATETIME returns no fractional seconds in Date!
%
% see also: DATENUM, DATEVEC
%

d      = datevec(t);

d(:,6) = round(d(:,6));  % Round seconds

dd = d(:,3);  % Original DayNumber

quot = [ 60 60 24 ]; %  [ ss-->mm  mm-->hh  hh-->dd ]
ind  = [ 6  5  4  ];

for ii = 1 : 3

    p = fix( d(:,ind(ii)) / quot(ii) );

 d(:,ind(ii)-0) = d(:,ind(ii)-0) - p * quot(ii);
 d(:,ind(ii)-1) = d(:,ind(ii)-1) + p;
  
end

% Check if DayNumber has changed

ii = find( d(:,3) > dd );

if isempty(ii)
   d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
   return
end

% New Date

[d(ii,1),d(ii,2),d(ii,3)] = datevec( datenum(d(ii,1),d(ii,2),d(ii,3)) );


d(:,1:3) = d(:,1:3) + all( d(:,1:3) == 0  , 2 ) * [ -1 12 31 ];
