classdef (Sealed, InferiorClasses = {?duration, ?calendarDuration ?matlab.graphics.axis.Axes}) datetime < matlab.mixin.internal.MatrixDisplay
%DATETIME Arrays to represent dates and times.
%   datetime arrays store values that represent points in time, including a date
%   and a time of day. Use the DATETIME constructor to create an array of datetimes
%   from strings, character vectors, or from vectors of date/time components. 
%   Use DATETIME('now'), DATETIME('today'), DATETIME('yesterday'), or DATETIME('tomorrow') 
%   to create scalar datetimes at or around the current moment.
%
%   You can subscript and manipulate datetime arrays just like ordinary numeric
%   arrays. Datetime arrays also support sorting and comparison, mathematical
%   calculations, as well as operations involving date and time components.
%
%   Each element of a datetime array represents one point in time. Use a
%   duration array to represent lengths of time in fixed-length time units.
%   Use a calendarDuration array to represent lengths of time in terms of
%   flexible-length calendar units.
%                       
%   A datetime array T has properties that store metadata, such as the display
%   format, and properties that allow you to access and modify the array's
%   values via its date/time components.  Access or assign to a property using
%   P = T.PropName or T.PropName = P, where PropName is one of the following:
%
%   DATETIME properties:
%       Format   - A character vector describing the format in which the array's values
%                  display.
%       TimeZone - A character vector representing the time zone in which the array's
%                  values are interpreted.
%       Year     - An array containing each element's year number.
%       Month    - An array containing each element's month number.
%       Day      - An array containing each element's day of month number.
%       Hour     - An array containing each element's hour.
%       Minute   - An array containing each element's minute.
%       Second   - An array containing each element's second, including a
%                  fractional part.
%
%   DATETIME methods and functions:
%     Creating arrays of datetimes:
%       datetime           - Create an array of datetimes.
%       isdatetime         - True for a array of datetimes.
%     Extract date and time components:
%       ymd                - Year, month, and day numbers of datetimes.
%       hms                - Hour, minute, and second numbers of datetimes.
%       year               - Year numbers of datetimes.
%       quarter            - Quarter numbers of datetimes.
%       month              - Month numbers or names of datetimes.
%       week               - Week numbers of datetimes.
%       day                - Day numbers or names of datetimes.
%       hour               - Hour numbers of datetimes.
%       minute             - Minute numbers of datetimes.
%       second             - Second numbers of datetimes.
%       timeofday          - Elapsed time since midnight for datetimes.
%       tzoffset           - Time zone offset of datetimes.
%       isdst              - True for datetimes occurring during Daylight Saving Time.
%       isweekend          - True for datetimes occurring on a weekend.
%     Calendar calculations with datetimes:
%       dateshift          - Shift datetimes or generate sequences according to a calendar rule.
%       between            - Difference between datetimes as calendar durations.
%       caldiff            - Successive differences between datetimes as calendar durations.
%     Mathematical calculations with datetimes:
%       plus               - Datetime addition.
%       minus              - Datetime subtraction.
%       diff               - Successive differences between datetimes as durations.
%       colon              - Create equally-spaced sequence of datetimes.
%       linspace           - Create equally-spaced sequence of datetimes.
%       mean               - Mean of datetimes.
%       median             - Median of datetimes.
%       mode               - Most frequent datetime value.
%       isnat              - True for datetimes that are Not-a-Time.
%       isinf              - True for datetimes that are +Inf or -Inf.
%       isfinite           - True for datetimes that are finite.
%     Comparisons between datetimes:
%       eq                 - Equality comparison for datetimes.
%       ne                 - Not-equality comparison for datetimes.
%       lt                 - Less than comparison for datetimes.
%       le                 - Less than or equal comparison for datetimes.
%       ge                 - Greater than or equal comparison for datetimes.
%       gt                 - Greater than comparison for datetimes.
%       isbetween          - Determine if datetimes are contained in an interval.
%       min                - Find minimum of datetimes.
%       max                - Find maximum of datetimes.
%       sort               - Sort datetimes.
%       sortrows           - Sort rows of a datetime array.
%       issorted           - True for sorted datetime vectors and matrices.
%     Set membership:
%       intersect          - Find datetimes common to two arrays.
%       ismember           - Find datetimes in one array that occur in another array.
%       setdiff            - Find datetimes that occur in one array but not in another.
%       setxor             - Find datetimes that occur in one or the other of two arrays, but not both.
%       unique             - Find unique datetimes in an array.
%       union              - Find datetimes that occur in either of two arrays.
%     Plotting:
%       plot               - Plot datetimes.
%     Conversion to other numeric representations:
%       exceltime          - Convert datetimes to Excel serial day numbers.
%       posixtime          - Convert datetimes to Posix time values.
%       juliandate         - Convert datetimes to Julian dates.
%       yyyymmdd           - Convert datetimes to YYYYMMDD numeric values.
%       datenum            - Convert datetimes to datenum values.
%       datevec            - Convert datetimes to date vectors.
%     Conversion to strings:
%       cellstr            - Convert datetimes to cell array of character vectors.
%       char               - Convert datetimes to character matrix.
%       datestr            - Convert datetimes to character vectors.
%       string             - Convert datetimes to strings.
%
%   Examples:
%
%      % Create a datetime array for the first 5 months in 2014.
%      t1 = datetime(2014,1:5,1)
%
%      % Add a random number of calendar days to each datetime. Extract
%      % the day component.
%      t1 = t1 + caldays(randi([0 15],1,5))
%      day = t1.Day
%
%      % Add a random amount of time to each datetime.
%      t1 = t1 + hours(rand(1,5))
%
%      % Shift each datetime to the end of its month.
%      t2 = dateshift(t1,'end','month')
%
%      % Find the time difference in hours/minutes/seconds between the two
%      % sets of datetimes.
%      d = t2 - t1
%
%      % Find the calendar time difference between the two sets of datetimes.
%      d2 = between(t1,t2)
%
%   See also DATETIME, DURATION.

%   Copyright 2014-2016 The MathWorks, Inc.
    
    properties(GetAccess='public', Dependent=true)
        %FORMAT Display format property for datetime arrays.
        %   The Format property specifies the format used to display the datetimes in
        %   the array. This property is a character vector constructed using the characters A-Z
        %   and a-z to represent date and time components of the datetimes. See the
        %   <a href="matlab:doc('datetime.Format')">datetime.Format property reference page</a> for the complete specification.
        %
        %   Changing the display format does not change the datetime values in the
        %   array, only their display.
        %
        %   The factory setting for the default value when you create a datetime array
        %   is locale-dependent. For information on how to change the
        %   default in the Preferences dialog box, see <a href="matlab:helpview(fullfile(docroot,'matlab/matlab_env/command-window-preferences.html'))">Set Command Window Preferences</a>.
        %   Datetime arrays whose time zone is set to 'UTCLeapSeconds' must
        %   use 'uuuu-MM-dd''T''HH:mm:ss[.SSS]Z', where you may specify
        %   from 0 to 9 fractional seconds digits.
        %
        %   See also DATETIME.
        Format
        
        %TIMEZONE Time zone property for datetime arrays.
        %   The TimeZone property array specifies the time zone used to interpret the
        %   datetimes in the array. Specify the time zone as:
        %
        %      - '' to create "unzoned" datetimes that do not belong to a specific
        %        time zone.
        %      - The name of a time zone region from the IANA Time Zone Database, e.g.
        %        'America/Los_Angeles'. The array obeys the time zone offset and Daylight
        %        Saving Time rules associated with that region.
        %      - An ISO 8601 character vector of the form +HH:MM or -HH:MM.
        %      - 'UTC' to create datetimes in Universal Coordinated Time.
        %      - 'UTCLeapSeconds' to create datetimes in Universal Coordinated Time that
        %        account for leap seconds.
        %
        %   The default value for TimeZone when you create a datetime array is ''.
        %   Datetime arrays with no time zone can not be compared or combined with
        %   arrays that have their TimeZone property set to a specific time zone.
        %
        %   Changing the TimeZone property of a datetime array from one time zone to
        %   another does not change the underlying points in time that the array's
        %   elements represent. Only the representation in terms of days, hours, etc.
        %   changes. Changing the TimeZone property from '' to a specific time zone puts
        %   the datetime values in that time zone without altering their Year, Month,
        %   Day, Hour, Minute, and Second proerties.
        %
        %   See also DATETIME, TIMEZONES.
        TimeZone
        
        %YEAR Datetime array year property.
        %   The Year property contains the year number of each datetime in the array.
        %   This property is the same size and shape as the datetime array.
        %
        %   Each year number is an integer value based on the proleptic Gregorian
        %   calendar. Years in the current era are positive, years in the previous era
        %   are zero or negative. For example, the year number of 1 BCE is 0.
        %
        %   If you set the Year property to a non-leap year for a datetime that occurs
        %   on a leap day (Feb 29th), the Day and Month properties change to Mar 1st.
        %
        %   See also DATETIME, MONTH, DAY, HOUR, MINUTE, SECOND
        Year

        %MONTH Datetime array month property.
        %   The Month property contains the month number of each datetime in the array.
        %   This property is the same size and shape as the datetime array.
        %
        %   Each month number is an integer value from 1 to 12, based on the proleptic
        %   Gregorian calendar. If you set a value outside that range, the Year property
        %   adjusts accordingly, and the Month property stays within that range. For
        %   example, month 0 corresponds to month 12 of the previous year.
        %
        %   If you change the Month property and the existing value of the Day property
        %   exceeds the length of the new month, the Month and Day property adjust
        %   accordingly. For example, if you set the Month property to 4 for a datetime
        %   that occurs on Jan 31, the resulting datetime adjusts to Jun 1.
        %
        %   See also DATETIME, YEAR, DAY, HOUR, MINUTE, SECOND
        Month

        %DAY Datetime array day property.
        %   The Day property contains the day number of each datetime in the array.
        %   This property is the same size and shape as the datetime array.
        %
        %   Each day number is an integer value from 1 to 28, 29, 30, or 31, depending
        %   on the month and year, and is based on the proleptic Gregorian calendar. If
        %   you set a value outside that range, the Month and Year properties adjust
        %   accordingly, and the Day property stays within that range. For example, day
        %   0 corresponds to the last day of the previous month.
        %
        %   See also DATETIME, YEAR, MONTH, HOUR, MINUTE, SECOND
        Day

        %HOUR Datetime array Hour property.
        %   The Hour property contains the hour number of each datetime in the array.
        %   This property is the same size and shape as the datetime array.
        %
        %   Each hour number is an integer value from 0 to 23. If you set a value
        %   outside that range, the Day, Month, and Year properties adjust accordingly,
        %   and the Hour property stays within that range. For example, hour -1
        %   corresponds to hour 23 of the previous day.
        %
        %   If a value you use to set the Hour property would create a non-existent
        %   datetime in the "spring ahead" gap of a daylight saving time shift, the Hour
        %   property adjusts to the next hour. If a value you use to set the Hour
        %   property would create an ambiguous datetime in the "fall back" overlap of a
        %   daylight saving time shift, the datetime adjusts to the second of the two
        %   times (i.e. in standard time) with that hour.
        %
        %   See also DATETIME, YEAR, MONTH, DAY, MINUTE, SECOND
        Hour

        %MINUTE Datetime array Minute property.
        %   The Minute property contains the minute number of each datetime in the
        %   array. This property is the same size and shape as the datetime array.
        %
        %   Each minute number is an integer value from 0 to 59. If you set a value
        %   outside that range, the Hour, Day, Month, and Year properties adjust
        %   accordingly, and the Minute property stays within that range. For example,
        %   minute -1 corresponds to minute 59 of the previous hour.
        %
        %   See also DATETIME, YEAR, MONTH, DAY, HOUR, SECOND
        Minute

        %SECOND Datetime array Second property.
        %   The Second property contains the second of each datetime in the array.
        %   This property is the same size and shape as the datetime array.
        %
        %   Each second value is a floating point value ordinarily ranging from 0 to
        %   strictly less than 60. If you set a value outside that range, the Minute,
        %   Hour, Day, Month, and Year properties adjust accordingly, and the Second
        %   property stays within that range. For example, second -1 corresponds to
        %   second 59 of the previous minute.
        %
        %   However, a datetime array whose TimeZone property is set to 'UTCLeapSeconds'
        %   has seconds ranging from 0 to strictly less than 61, with values from 60 to 
        %   61 for datetimes that are during a leap second occurrence.
        %
        %   A datetime array represents points in time to an accuracy of at least 1 ns.
        %
        %   See also DATETIME, YEAR, MONTH, DAY, HOUR, MINUTE
        Second
    end
    properties(GetAccess='public', Constant)
    %SYSTEMTIMEZONE System time zone setting.
    %   The SystemTimeZone property contains the time zone that the system is set to.
    %
    %   See also TIMEZONE.
        SystemTimeZone = matlab.internal.datetime.getDefaults('SystemTimeZone');
    end
    properties(GetAccess='public', Hidden, Constant)
        % These properties are for internal use only and will change in a
        % future release.  Do not use these properties.
        UTCZoneID = 'UTC';
        UTCLeapSecsZoneID = 'UTCLeapSeconds';
        ISO8601Format = 'uuuu-MM-dd''T''HH:mm:ss.SSS''Z''';
        epochDN = 719529; % 1-Jan-1970 00:00:00
        
        % The MonthsOfYear property contains the month names localized in
        % the system locale. The property is a scalar struct with 'Short'
        % and 'Long' fields, each containing a cell array of character vectors for
        % the short and long month names, respectively.
        MonthsOfYear = struct('Short',{matlab.internal.datetime.getMonthNames('short', getDatetimeSettings('locale'))}, ...
                              'Long',{matlab.internal.datetime.getMonthNames('long', getDatetimeSettings('locale'))});
                              
        % The DayNames property contains the day names localized in the system locale.
        % The property is a scalar struct with 'Short' and 'Long' fields, each containing
        % a cell array of character vectors for the short and long day names, respectively.
        DaysOfWeek = struct('Short',{matlab.internal.datetime.getDayNames('short', getDatetimeSettings('locale'))}, ...
                            'Long',{matlab.internal.datetime.getDayNames('long', getDatetimeSettings('locale'))});
    end
    
    properties(GetAccess='protected', SetAccess='protected')
        % Count seconds (including fractional) from epoch as a double-double
        % stored in complex
        data
        
        % Display format.  When this field is empty the global setting will 
        % be used
        fmt = '';
        
        % A time zone ID, e.g. America/Los_Angeles. Empty character vector means the array
        % contains no time zone information, i.e., it's not a fully specified
        % value, and cannot be mixed with "real-world" datetimes
        tz
        
        % Determines whether or not we are displaying the time when
        % displaying datetime arrays
        isDateOnly = false;
    end
    
    methods(Access = 'public')
        function this = datetime(inData,varargin)
        %DATETIME Create an array of datetimes.
        %   D = DATETIME('now') returns the current date and time. D is a scalar
        %   datetime. A datetime is a value that represents a point in time. D is
        %   "unzoned", i.e., does not belong to a particular time zone. D = DATETIME
        %   with no inputs is equivalent to D = DATETIME('now').
        %
        %   D = DATETIME('today') returns the current date. D is an unzoned scalar
        %   datetime with the time portion set to 00:00:00.
        %
        %   D = DATETIME('tomorrow') returns the date of the following day. D is an
        %   unzoned scalar datetime with the time portion set to 00:00:00.
        %
        %   D = DATETIME('yesterday') returns the date of the previous day. D is an
        %   unzoned scalar datetime with the time portion set to 00:00:00.
        %
        %   D = DATETIME(DS,'InputFormat',INFMT) creates an array of datetimes from the
        %   character vectors in the cell array DS, or the strings in the string array DS. 
        %   DS contains date/time strings or character vectors in the format specified by 
        %   the format string INFMT. D also displays using INFMT. All
        %   strings in DS must have the same format.
        %
        %   INFMT is a datetime format string constructed using the characters A-Z and
        %   a-z to represent date and time components of the strings. See the description
        %   of the <a href="matlab:doc('datetime.Format')">Format property</a> for details. If you do not provide INFMT, DATETIME
        %   tries to determine the format automatically by first trying the value of the
        %   'Format' parameter, followed by the default display format, followed by some
        %   other common formats.
        %
        %   For best performance, and to avoid ambiguities between MM/dd and dd/MM formats,
        %   always specify INFMT.
        %
        %   If INFMT does not include a date portion, DATETIME assumes the current day.
        %   If INFMT does not include a time portion, DATETIME assumes midnight. DATETIME
        %   creates unzoned datetimes and ignores any time zone information in DS when
        %   INFMT includes a time zone specifier.
        %
        %   D = DATETIME(DS,'InputFormat',INFMT,'Locale',LOCALE) specifies the locale that
        %   DATETIME uses to interpret the strings or character vectors in DS. LOCALE must be a 
        %   character vector in the form xx_YY, where xx is a lowercase ISO 639-1 two-letter 
        %   language code and YY is an uppercase ISO 3166-1 alpha-2 country code, for example 
        %   ja_JP. LOCALE can also be the character vector 'system' to use the system locale setting. 
        %   All strings or character vectors in DS must have the same locale.
        %
        %   D = DATETIME(DS,'InputFormat',INFMT,'PivotYear',PIVOT,...) specifies the pivot
        %   year that DATETIME uses to interpret the strings or character vectors in DS when FMT 
        %   contains the two-digit year specifier 'yy'. The default pivot is YEAR(DATETIME('NOW')) - 50.
        %   PIVOT only affects the year specifier 'yy' in FMT.
        %
        %   D = DATETIME(DV) creates a column vector of datetimes from a numeric matrix
        %   DV in the form [YEAR MONTH DAY HOUR MINUTE SECOND]. The first five columns
        %   must contain integer values, while the last may contain fractional seconds.
        %   DV may also be in the form [YEAR MONTH DAY], in which case the hours,
        %   minutes, and seconds are taken to be zero.
        %
        %   D = DATETIME(Y,MO,D,H,MI,S) or DATETIME(Y,MO,D) creates an array of
        %   datetimes from separate arrays. The arrays must be the same size, or any can
        %   be a scalar. Y, MO, D, H, and MI must contain integer values, while S may
        %   contain fractional seconds.
        %
        %   D = DATETIME(Y,MO,D,H,MI,S,MS) creates an array of datetimes from separate
        %   arrays, including milliseconds. The arrays must be the same size, or any can
        %   be a scalar. Y, MO, D, H, MI, and S must contain integer values, while MS
        %   may contain fractional milliseconds.
        %
        %   D = DATETIME(X,'ConvertFrom',TYPE) converts the numeric values in X to a
        %   DATETIME array D. D is the same size as X. TYPE specifies the type of values
        %   contained in X, and is one of the following character vectors:
        %
        %      'datenum'             The number of days since 0-Jan-0000 (Gregorian)
        %      'posixtime'           The number of seconds since 1-Jan-1970 00:00:00 UTC
        %      'excel'               The number of days since 0-Jan-1900
        %      'excel1904'           The number of days since 0-Jan-1904
        %      'juliandate'          The number of days since noon UTC 24-Nov-4714 BCE (Gregorian)
        %      'modifiedjuliandate'  The number of days since midnight UTC 17-Nov-1858
        %      'yyyymmdd'            A YYYYMMDD numeric value, e.g. 20140716
        %
        %   D = DATETIME(X,'ConvertFrom','epochtime','Epoch',EPOCH) converts the numeric
        %   values in X to a DATETIME array D. X contains the number of seconds since
        %   the epoch time. EPOCH is a scalar DATETIME or a date/time character vector, representing
        %   the epoch time.
        %
        %   D = DATETIME(...,'Format',FMT) creates D with the specified display format.
        %   FMT is a datetime format character vector. See the description of the <a href="matlab:doc('datetime.Format')">Format property</a>
        %   for details. The factory setting for the default is locale-dependent. For
        %   information on setting the default in the Preferences dialog box, see 
        %   <a href="matlab:helpview(fullfile(docroot,'matlab/matlab_env/command-window-preferences.html'))">Set Command Window Preferences</a>.  
        %   When creating D from strings or character vectors, specify
        %   FMT as 'preserveinput' to use the the 'InputFormat' parameter
        %   (or the format that was determined automatically if
        %   'InputFormat' was not given). Specify FMT as 'default' to use
        %   the default display format.
        %
        %   D = DATETIME(...,'TimeZone',TZ,...) specifies the time zone used to interpret
        %   the input data, and that the datetimes in D are in. TZ is the name of a time
        %   zone region, e.g. 'America/Los_Angeles', or an ISO 8601 character vector of the form
        %   +HH:MM or -HH:MM. Specify 'UTC' to create datetimes in Universal Coordinated
        %   Time. Specify the empty character vector '' to create "unzoned" datetimes that do not
        %   belong to a specific time zone. Specify 'local' to create datetimes in the
        %   system time zone. See the description of the <a href="matlab:doc('datetime.TimeZone')">TimeZone property</a> for more
        %   details. When the input data are strings or character vectors that include time zone offsets such
        %   as 'EST' or 'CEST', DATETIME converts all datetimes to the time zone specified
        %   by TZ.
        %
        %   Examples:
        %
        %      % Create a scalar datetime representing the current date and time.
        %      t = datetime
        %
        %      % Create a datatime array for the first 5 days in July, 2014.
        %      t = datetime(2014,7,1:5)
        %
        %      % Create a datetime array for the 5 days following 28 July, 2014.
        %      t = datetime(2014,7,28+(1:5))
        %
        %      % Create a datetime array from character vectors. Have the array use the default
        %      % format for display, then have it use the original format from the character vectors.
        %      s = {'2014-07-28' '2014-07-29' '2014-07-30'}
        %      t1 = datetime(s,'InputFormat','yyyy-MM-dd')
        %      t2 = datetime(s,'Format','yyyy-MM-dd')
        %
        %   See also DATETIME, DURATION

            import matlab.internal.datetime.createFromString
            import matlab.internal.table.parseArgs
            
            if nargin == 0 % same as datetime('now')
                thisData = currentTime();
                fmt = '';
                this.isDateOnly = false;
                tz = '';
            
            else
                haveStrings = false;
                haveNumeric = false;
                haveNamedInstant = false;
            
                if ischar(inData)
                    [haveNamedInstant,yesterdayTodayTomorrow] = ismember(lower(inData),{'now' 'yesterday' 'today' 'tomorrow'});
                    if ~haveNamedInstant
                        haveStrings = true;
                        inData = cellstr(inData);
                    end
                elseif isnumeric(inData)
                    haveNumeric = true;
                elseif isa(inData,'string')
                    haveStrings = true;                                        
                elseif matlab.internal.datatypes.isCharStrings(inData)
                    haveStrings = true;
                end
                
                % Find how many numeric inputs args: count up until the first non-numeric.
                numNumericArgs = haveNumeric + sum(cumprod(cellfun(@isnumeric,varargin)));

                dfltFmt = '';
                this.isDateOnly = false;
                tryFmts = {};
                dfltPivot = matlab.internal.datetime.getDefaults('pivotyear');
                preserveInputFmt = false;
                if isempty(varargin)
                    % Default format and the local time zone and locale.
                    inputFmt = dfltFmt;
                    fmt = dfltFmt; %#ok<*PROP>
                    tz = '';
                    locale = '';
                    pivot = dfltPivot;
                    epoch = [];
                    supplied = struct('ConvertFrom',false, 'InputFormat',false, 'Format',false, 'TimeZone',false, 'Locale',false, 'PivotYear',false, 'Epoch',false);
                    isUTCLeapSecs = false;
                else
                    % Process explicit parameter name/value pairs.
                    pnames = {'ConvertFrom' 'InputFormat' 'Format' 'TimeZone' 'Locale'      'PivotYear' 'Epoch'};
                    dflts =  {''            ''             dfltFmt ''          ''            dfltPivot   0     };
                    args = varargin(max(numNumericArgs,1):end);
                    [convertFrom,inputFmt,fmt,tz,locale,pivot,epoch,supplied] = parseArgs(pnames, dflts, args{:});

                    % Canonicalize the TZ name. This make US/Eastern ->
                    % America/New_York, but it will also make EST -> Etc/GMT+5,
                    % because EST is an offset, not a time zone.
                    if supplied.TimeZone, tz = verifyTimeZone(tz); end
                    if supplied.Locale, locale = matlab.internal.datetime.verifyLocale(locale); end
                    if supplied.PivotYear, verifyPivot(pivot); end
                    isUTCLeapSecs = strcmp(tz,datetime.UTCLeapSecsZoneID);
                    if isUTCLeapSecs
                        if supplied.InputFormat
                            if ~strcmp(inputFmt,datetime.ISO8601Format)
                                error(message('MATLAB:datetime:InvalidUTCLeapSecsFormatString',datetime.ISO8601Format));
                            end
                        else
                            inputFmt = datetime.ISO8601Format;
                        end
                        if supplied.Format
                            if ~strcmp(fmt,datetime.ISO8601Format)
                                error(message('MATLAB:datetime:InvalidUTCLeapSecsFormatString',datetime.ISO8601Format));
                            end
                        else
                            fmt = datetime.ISO8601Format;
                        end
                    else
                        if supplied.Format
                            if strcmpi(fmt,'preserveinput')
                                preserveInputFmt = true;
                                supplied.Format = false;
                            else
                                checkConflicts = haveStrings && ~supplied.InputFormat; % only when also used to read strings
                                [fmt,this.isDateOnly] = verifyFormat(fmt,tz,checkConflicts,true);
                            end
                        end
                    end
                end
                
                if isa(inData,'datetime') % construct from an array of datetimes
                    % Take values from the input array rather than from the defaults.
                    if ~supplied.Format, fmt = inData.fmt; end
                    if ~supplied.TimeZone, tz = inData.tz; end

                    % Adjust for a new time zone.
                    thisData = timeZoneAdjustment(inData.data,inData.tz,tz);
                    
                elseif haveNumeric
                    if supplied.ConvertFrom % datetime(x,'ConvertFrom',type,...)
                        if numNumericArgs > 1
                            error(message('MATLAB:datetime:WrongNumInputsConversion'));
                        end
                        try %#ok<ALIGN>
                            thisData = datetime.convertFrom(inData,convertFrom,tz,epoch);
                        catch ME, throw(ME), end
                    else
                        if numNumericArgs == 1 % datetime([y,mo,d],...) or datetime([y,mo,d,h,mi,s],...)
                            ncols = size(inData,2);
                            if ~ismatrix(inData) || ((ncols ~= 3) && (ncols ~= 6))
                                error(message('MATLAB:datetime:InvalidNumericData'));
                            end
                            inData = num2cell(full(double(inData)),1); % split into separate vectors.
                        else % datetime(y,mo,d,...), datetime(y,mo,d,h,mi,s,...), or datetime(y,mo,d,h,mi,s,ms,...)
                            inData = expandNumericInputs([{inData} varargin(1:numNumericArgs-1)]);
                        end
                        if ~supplied.Format && (numel(inData) == 3) && ~isUTCLeapSecs
                            this.isDateOnly = true;
                        end
                    
                        try %#ok<ALIGN>
                            thisData = matlab.internal.datetime.createFromDateVec(inData,tz); % or datevec + millis
                        catch ME, throw(ME), end
                    end
                    
                elseif haveNamedInstant
                    % Get the system clock. If the requested result is zoned, the system time zone's
                    % offset is removed to translate the value to UTC.
                    thisData = currentTime(tz);
                    if yesterdayTodayTomorrow > 1
                        % Floor to get today, then subtract or add a day for yesterday or tomorrow.
                        if ~supplied.Format, this.isDateOnly = true; end
                        ucal = getDateFieldsStructure;
                        thisData = matlab.internal.datetime.datetimeFloor(thisData,ucal.DAY_OF_MONTH,tz);
                        thisData = matlab.internal.datetime.addToDateField(thisData,yesterdayTodayTomorrow-3,ucal.DAY_OF_MONTH,tz);
                    end
                
                elseif haveStrings
                    % Construct from a cell array of date strings of any shape. Error if none
                    % of the strings can be parsed. Returning all NaT would give an indication
                    % that something went wrong, but offer s no help on what to do. For example,
                    % error if one string is given and it can't be parsed.                   
                    if isUTCLeapSecs
                        try %#ok<ALIGN>
                            thisData = createFromString(inData,inputFmt,1,tz,locale,pivot);
                        catch ME, throw(ME), end
                    else
                        try
                            if supplied.InputFormat
                                verifyFormat(inputFmt,tz,true);
                                thisData = createFromString(inData,inputFmt,1,tz,locale,pivot);
                            else
                                if supplied.Format, tryFmts = {fmt}; end
                                [thisData,inputFmt] = guessFormat(inData,tryFmts,1,tz,locale,pivot);
                            end
                        catch ME
                            if strcmp(ME.identifier,'MATLAB:datetime:ParseErrs')
                                handleParseErrors(inData,supplied,fmt,inputFmt,locale);
                            else
                                throw(ME)
                            end
                        end
                        if ~supplied.Format
                            if preserveInputFmt
                                fmt = inputFmt;
                            elseif isDateOnlyFormat(inputFmt)
                                this.isDateOnly = true;
                            end
                        end
                    end
                    
                else
                    error(message('MATLAB:datetime:InvalidData'));
                end
            end
            
            this.data = thisData;
            this.fmt = fmt;
            this.tz = tz;
        end
        
        %% Conversions to numeric types
        % No direct conversions, need to subtract a time origin
        
        %% Conversions to string types
        function s = char(this,format,locale)
            %CHAR Convert datetimes to character matrix.
            %   C = CHAR(T) returns a character matrix representing the datetimes in T.
            %
            %   C = CHAR(T,FMT) uses the specified datetime format. FMT is a character vector
            %   constructed using the characters A-Z and a-z to represent date and time
            %   components of the datetimes. See the description of the <a href="matlab:doc('datetime.Format')">Format property</a>
            %   for details.
            %
            %   C = CHAR(T,FMT,LOCALE) specifies the locale (in particular, the language)
            %   used to create the character vectors. LOCALE must be a character vector in the form
            %   xx_YY, where xx is a lowercase ISO 639-1 two-letter language code and YY is
            %   an uppercase ISO 3166-1 alpha-2 country code, for example 'ja_JP'. LOCALE
            %   can also be the character vector 'system' to use the system locale setting. 
            %
            %   See also CELLSTR, STRING, DATETIME.
            if nargin < 2 || isequal(format,[])
                format = getDisplayFormat(this);
            else
                format = verifyFormat(format,this.tz,false,true);
            end
            if nargin < 3 || isequal(locale,[])
                s = char(matlab.internal.datetime.formatAsString(this.data,format,this.tz,false,getDatetimeSettings('locale')));
            else
                s = char(matlab.internal.datetime.formatAsString(this.data,format,this.tz,false,matlab.internal.datetime.verifyLocale(locale)));
            end
        end
        
        function c = cellstr(this,format,locale)
            %CELLSTR Convert datetimes to cell array of character vectors.
            %   C = CELLSTR(T) returns a cell array of character vectors representing the datetimes in T.
            %
            %   C = CELLSTR(T,FMT) uses the specified datetime format. FMT is a character vector
            %   constructed using the characters A-Z and a-z to represent date and time
            %   components of the datetimes. See the description of the <a href="matlab:doc('datetime.Format')">Format property</a>
            %   for details.
            %
            %   C = CELLSTR(T,FMT,LOCALE) specifies the locale (in particular, the language)
            %   used to create the character vectors. LOCALE must be a character vector in the form xx_YY, where
            %   xx is a lowercase ISO 639-1 two-letter language code and YY is an uppercase
            %   ISO 3166-1 alpha-2 country code, for example 'ja_JP'. LOCALE can also be the
            %   character vector 'system' to use the system locale setting. 
            %
            %   See also CHAR, STRING, DATETIME.
            if nargin < 2 || isequal(format,[])
                format = getDisplayFormat(this);
            else
                format = verifyFormat(format,this.tz,false,true);
            end
            if nargin < 3 || isequal(locale,[])
                c = matlab.internal.datetime.formatAsString(this.data,format,this.tz,false,getDatetimeSettings('locale'));
            else
                c = matlab.internal.datetime.formatAsString(this.data,format,this.tz,false,matlab.internal.datetime.verifyLocale(locale));
            end
        end
        
        function s = string(this,format,locale)
            %STRING Convert datetimes to strings.
            %   S = STRING(T) returns a string array representing the datetimes in T.
            %
            %   S = STRING(T,FMT) uses the specified datetime format. FMT is a char vector
            %   constructed using the characters A-Z and a-z to represent date and time
            %   components of the datetimes, for example 'dd-MMM-yyyy'. See the description
            %   of the <a href="matlab:doc('datetime.Format')">Format property</a> for details.
            %
            %   S = STRING(T,FMT,LOCALE) specifies the locale (in particular, the language)
            %   used to create the strings. LOCALE must be a char vector in the form xx_YY,
            %   where xx is a lowercase ISO 639-1 two-letter language code and YY is an
            %   uppercase ISO 3166-1 alpha-2 country code, for example 'ja_JP'. LOCALE can
            %   also be the character vector 'system' to use the system locale setting. 
            %
            %   See also CHAR, CELLSTR, DATETIME.
            if nargin < 2 || isequal(format,[])
                format = getDisplayFormat(this);
            else
                format = verifyFormat(format,this.tz,false,true);
            end
            if nargin < 3 || isequal(locale,[])
                s = matlab.internal.datetime.formatAsString(this.data,format,this.tz,true,getDatetimeSettings('locale'));
            else
                s = matlab.internal.datetime.formatAsString(this.data,format,this.tz,true,matlab.internal.datetime.verifyLocale(locale));
            end
        end
        
        %% Conversions to the legacy types
        function dn = datenum(this)
            %DATENUM Convert datetimes to serial date numbers.
            %   DN = DATENUM(T) converts the datetimes in the array T to the equivalent
            %   MATLAB serial date number, i.e. the number of days and fractional days
            %   since 0-Jan-0000 00:00:00. D is a double array.
            %
            %   DATENUM treats T as if it were unzoned. Elements of T that occur during
            %   Daylight Saving Time are not adjusted. For example, datetimes that are 1
            %   hour apart in real time, but that span a "spring forward" DST shift have
            %   date numbers that differ by 2/24, not 1/24.
            %
            %   See also EXCELTIME, POSIXTIME, JULIANDAY, DATEVEC, DATESTR, DATETIME.
            
            % Convert to unzoned, no leap seconds.
            thisData = timeZoneAdjustment(this.data,this.tz,'');
            
            % Get the day number (including fractional days) since 0-Jan-0000 00:00:00.
            millisPerDay = 86400*1000;
            datenumOffset = datetime.epochDN*millisPerDay;
            dn = (real(thisData) + datenumOffset) / millisPerDay; % round trip exact up to ms
        end
        
        function s = datestr(this,varargin)
            %DATESTR Convert datetimes to character matrix.
            %   C = DATESTR(T) returns a character matrix representing the
            %   datetimes in T.
            %
            %   S = datestr(T,F) converts the datetimes in T using the format F.
            %   NOTE: valid values for F are different than those accepted for the
            %   Format property of T. Type "help datestr" for a description of
            %   valid values for F.
            %
            %   See also CELLSTR, CHAR, DATENUM, DATEVEC, DATETIME.
            s = datestr(datenum(this),varargin{:});
        end
        
        function [y,mo,d,h,m,s] = datevec(this)
            %DATEVEC Convert datetimes to date vectors.
            %   DV = DATEVEC(T) splits the datetime array T into separate column vectors for
            %   years, months, days, hours, minutes, and seconds, and returns one numeric
            %   matrix.
            %
            %   [Y,MO,D,H,MI,S] = DATEVEC(T) returns the components of T as individual
            %   variables.
            %
            %   See also YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, YMD, HMS,
            %            TIMEOFDAY, DATENUM, DATESTR, DATETIME.
            
            % This preserves the zoned component values.
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.EXTENDED_YEAR ucal.MONTH ucal.DAY_OF_MONTH ucal.HOUR_OF_DAY ucal.MINUTE ucal.SECOND];
            [y,mo,d,h,m,s] = matlab.internal.datetime.getDateFields(this.data,fieldIDs,this.tz);
            if nargout <= 1
                y = [y(:),mo(:),d(:),h(:),m(:),s(:)];
            end
        end
        
        %% Conversions to other time types
        function jd = juliandate(this,kind)
            %JULIANDATE Convert datetimes to Julian dates.
            %   JD = JULIANDATE(T) converts the datetimes in the array T to the equivalent
            %   Julian dates, i.e. the number of days and fractional days since noon on Nov
            %   24, 4714 BCE in the proleptic Gregorian calendar (Jan 1, 4713 BCE in the
            %   proleptic Julian calendar). If T has a time zone, JULIANDATE computes JD
            %   with respect to UTC. JULIANDATE ignores leap seconds unless T's time zone
            %   is 'UTCLeapSeconds'. D is a double array.
            %
            %   MJD = JULIANDATE(T,'modifiedjuliandate') converts the datetimes in the array
            %   T to the equivalent modified Julian date, i.e. the number of days and
            %   fractional days since Nov 17, 1858 00:00:00. JULIANDATE(T,'juliandate')
            %   returns Julian dates, the default behavior.
            %
            %   The modified Julian date is equal to the Julian date minus 2400000.5.
            %
            %   See also EXCELTIME, POSIXTIME, YYYYMMDD, DATENUM, DATETIME.
            if nargin == 1
                kind = 1;
            else
                kind = getChoice(kind,{'juliandate' 'jd' 'modifiedjuliandate' 'mjd'},[1 1 2 2]);
            end
            % This does not try to account for JD(UTC) vs. JD(UT1) vs. JD(TT), it simply
            % accepts the date/time components as is. The one exception is that for a
            % datetime set to 'UTCLeapSeconds', the fractional part of the JD on a day
            % with a leap second is normalized by 86401, not 86400.
            ucal = getDateFieldsStructure;
            jd = matlab.internal.datetime.getDateFields(this.data,ucal.MODIFIED_JULIAN_DATE,this.tz); % returns MJD
            if (kind == 1)
                JDtoMJDoffset = 2400000.5; % MJD = JD - 2400000.5
                jd = jd + JDtoMJDoffset;
            end
        end
        
        function p = posixtime(this)
            %POSIXTIME Convert datetimes to Posix times.
            %   P = POSIXTIME(T) converts the datetimes in the array T to the equivalent
            %   Posix time, i.e. the number of seconds (including fractional) that have
            %   elapsed since 00:00:00 1-Jan-1970 UTC, ignoring leap seconds. D is a double
            %   array.
            %
            %   See also EXCELTIME, JULIANDATE, YYYYMMDD, DATENUM, DATETIME.
            
            millisPerSec = 1000;
            thisData = this.data;
            if strcmp(this.tz,datetime.UTCLeapSecsZoneID)
                % POSIX time doesn't count leap seconds, remove them.
                [thisData,isLeapSec] = removeLeapSeconds(thisData);
                % The POSIX time during a leap second is defined equal to the corresponding
                % time in 0th second of the next minute, so shift leap seconds ahead.
                thisData = matlab.internal.datetime.datetimeAdd(thisData,millisPerSec*isLeapSec);
            end
            p = real(thisData) / millisPerSec; % ms -> s
        end
        
        function e = exceltime(this,timeSystem)
            %EXCELTIME Convert datetimes to Excel serial date numbers.
            %   E = EXCELTIME(T) converts the datetimes in the array T to the equivalent
            %   Excel serial date numbers, i.e. the number of days and fractional days since
            %   0-Jan-1900 00:00:00, ignoring time zone and leap seconds. D is a double
            %   array.
            %
            %   Excel serial date numbers treat 1900 as a leap year, therefore dates after
            %   28-Feb-1900 are off by 1 day relative to MATLAB serial date numbers, and
            %   there is a discontinuity of 1 day between 28-Feb-1900 and 1-Mar-1900.
            %
            %   E = EXCELTIME(T,'1904') converts the datetimes in the array T to the
            %   equivalent "1904-based" Excel serial date numbers, i.e. the number of days
            %   and fractional days since 1-Jan-1904 00:00:00, ignoring time zone.
            %   EXCELTIME(T,'1900') returns "1900-based" Excel serial date numbers, the
            %   default behavior.
            %
            %   EXCELTIME(T,'1904') is equal to EXCELTIME(T,'1900') - 1462.
            %
            %   NOTE: Excel serial date numbers are not defined prior to their epoch, i.e.
            %   prior to 0-Jan-1900 or 1-Jan-1904.
            %
            %   See also POSIXTIME, JULIANDATE, YYYYMMDD, DATENUM, DATETIME.
            thisData = timeZoneAdjustment(this.data,this.tz,'');
            millisPerDay = 86400*1000;
            excelOffset1900 = 25568 * millisPerDay;
            e = (real(thisData) + excelOffset1900) / millisPerDay; % consistent with datenum
            e = e + (e >= 60); % Correction for Excel's 1900 leap year bug
            if nargin > 1
                if strcmp(timeSystem,'1904') || isequal(timeSystem,1904)
                    e = e - 1462; % "1904" epoch is 0-Jan-1904
                elseif strcmp(timeSystem,'1900') || isequal(timeSystem,1900)
                    % OK
                else
                    error(message('MATLAB:datetime:exceltime:InvalidTimeSystem'));
                end
            end
        end
        
        function pd = yyyymmdd(this)
            %YYYYMMDD Convert MATLAB datetimes to YYYYMMDD numeric values.
            %   D = YYYYMMDD(T) returns a double array containing integers whose digits
            %   represent the datetime values in T. For example, the date July 16, 2014
            %   is converted to the integer 20140716. The conversion is performed as
            %   D = 10000*YEAR(T) + 100*MONTH(T) + DAY(T).
            %
            %   See also EXCELTIME, POSIXTIME, JULIANDATE, DATENUM, DATETIME.
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.EXTENDED_YEAR ucal.MONTH ucal.DAY_OF_MONTH];
            [y,mo,d] = matlab.internal.datetime.getDateFields(this.data,fieldIDs,this.tz);
            pd = y*10000 + mo*100 + d;
            
            % Preserve Infs. These have become NaNs from NaNs in month/day
            % components. Use year, which will be the appropriate non-finite.
            nonfinites = ~isfinite(pd);
            pd(nonfinites) = y(nonfinites);
        end
        
        
        %% Date/time component methods
        % These return datetime components as integer values (and sometimes
        % names), except for seconds, which returns non-integer values
        
        function [y,m,d] = ymd(this)
            %YMD Year, month, and day numbers of datetimes.
            %   [Y,M,D] = YMD(T) returns the year, month, and day numbers of the
            %   datetimes in T. Y, M, and D are numeric arrays the same size as
            %   T, containing integer values.
            %
            %   See also HMS, YEAR, QUARTER, MONTH, WEEK, DAY.
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.EXTENDED_YEAR ucal.MONTH ucal.DAY_OF_MONTH];
            [y,m,d] = matlab.internal.datetime.getDateFields(this.data,fieldIDs,this.tz);
        end
        
        function [h,m,s] = hms(this)
            %HMS Hour, minute, and second numbers of datetimes.
            %   [H,M,S] = HMS(T) returns the hour and minute numbers (as integer values) and
            %   the second values (including fractional part) of the datetimes in T. H, M,
            %   and S are numeric arrays the same size as T.
            %
            %   See also YMD, HOUR, MINUTE, SECOND.
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.HOUR_OF_DAY ucal.MINUTE ucal.SECOND];
            [h,m,s] = matlab.internal.datetime.getDateFields(this.data,fieldIDs,this.tz);
        end
        
        function y = year(this,kind)
            %YEAR Year numbers of datetimes.
            %   Y = YEAR(T) returns the year numbers of the datetimes in T. Y is an array
            %   the same size as T containing integer values.
            %
            %   Y = YEAR(T,KIND) returns the kind of year numbers specified by KIND. KIND
            %   is one of the character vectors
            %
            %      iso            - The ISO year number, which includes a year zero, and has
            %                       negative values for years BCE. This is the default.
            %      gregorian      - The Gregorian year number, which does not include a year
            %                       zero, and has positive values for years BCE.
            %
            % See also QUARTER, MONTH, WEEK, DAY, YMD.
            if nargin == 1
                kind = 1;
            else
                kind = getChoice(kind,{'ISO' 'Gregorian'},[1 2]);
            end
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.EXTENDED_YEAR ucal.YEAR];
            y = matlab.internal.datetime.getDateFields(this.data,fieldIDs(kind),this.tz);
        end
        
        function q = quarter(this)
            %QUARTER Quarter numbers of datetimes.
            %   Q = QUARTER(T) returns the quarter numbers of the datetimes in T. Q is an
            %   array the same size as T containing integer values from 1 to 4.
            %
            % See also YEAR, MONTH, WEEK, DAY, YMD.
            ucal = getDateFieldsStructure;
            q = matlab.internal.datetime.getDateFields(this.data,ucal.QUARTER,this.tz);
        end
        
        function m = month(this,kind)
            %MONTH Month numbers or names of datetimes.
            %   M = MONTH(T) returns the month numbers of the datetimes in T. M is an array
            %   the same size as T containing integer values from 1 to 12.
            %
            %   M = MONTH(T,KIND) returns the kind of month values specified by KIND. KIND
            %   is one of the character vectors
            %
            %      monthofyear - The month year number number. This is the default 
            %      name        - D is a cell array of character vectors the same size as T containing the
            %                    full month names. MONTH returns the empty character vector for NaT datetimes.
            %      shortname   - D is a cell array of character vectors the same size as T containing the
            %                    month name abbreviations. MONTH returns the empty character vector for NaT
            %                    datetimes.
            %
            % See also YEAR, QUARTER, WEEK, DAY, YMD.
            import matlab.internal.datetime.getMonthNames
            
            ucal = getDateFieldsStructure;
            fieldIDs = ucal.MONTH;
            m = matlab.internal.datetime.getDateFields(this.data,fieldIDs,this.tz);
            if nargin > 1
                kind = getChoice(kind,{'MonthOfYear' 'MoY' 'Name' 'LongName' 'ShortName'},[1 1 2 2 3]);
                if kind > 1
                    if kind == 2
                        names = matlab.internal.datetime.getMonthNames('long', getDatetimeSettings('locale'));
                    else % kind == 3
                        names = matlab.internal.datetime.getMonthNames('short', getDatetimeSettings('locale'));
                    end
                    names{end+1} = ''; % return empty character vector for NaT
                    m(isnan(m)) = length(names);
                    m = reshape(names(m),size(this));
                end
            end
        end
        
        function w = week(this,kind)
            %WEEK Week numbers of datetimes.
            %   W = WEEK(T) returns the week of year numbers of the datetimes in T. W is an
            %   array the same size as T containing integer values from 1 to 53.
            %
            %   W = WEEK(T,KIND) returns the kind of week numbers specified by KIND. KIND
            %   is one of the character vectors
            %
            %      weekofyear  - The week of year number. Jan 1st is defined to be in week 1 of its
            %                    year, even if fewer than 4 days of that week fall in the same year.
            %                    This is the default.
            %      weekofmonth - The week of month number, from 1 to 5. The 1st of the month is
            %                    defined to be in week 1 of its month, even if fewer than 4 days of
            %                    that week fall in the same month.
            %
            % See also YEAR, QUARTER, MONTH, DAY, YMD.
            if nargin == 1
                kind = 1;
            else
                kind = getChoice(kind,{'WeekOfYear' 'WoY' 'WeekOfMonth' 'WoM'},[1 1 2 2]);
            end
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.WEEK_OF_YEAR ucal.WEEK_OF_MONTH];
            w = matlab.internal.datetime.getDateFields(this.data,fieldIDs(kind),this.tz);
        end
        
        function d = day(this,kind) % DoM, DoY, DoW, ShortName, LongName
            %DAY Day numbers or names of datetimes.
            %   D = DAY(T) returns the day of month numbers of the datetimes in T. D is an
            %   array the same size as T containing integer values from 1 to 28, 29, 30, or
            %   31, depending on the month and year.
            %
            %   D = DAY(T,KIND) returns the kind of day values specified by KIND. KIND
            %   is one of the character vectors
            %
            %      dayofmonth - The day of month number number. This is the default 
            %      dayofweek  - The day of week number, from 1 to 7. Sunday is defined to be
            %                   day 1 of the week.
            %      dayofyear  - The day of year number, from 1 to 365 or 366, depending on the
            %                   year. This is sometimes incorrectly referred to as the "Julian
            %                   day".
            %      name       - D is a cell array of character vectors the same size as T containing the
            %                   full day names. DAY returns the empty character vector for NaT datetimes.
            %      shortname  - D is a cell array of character vectors the same size as T containing the
            %                   day name abbreviations. DAY returns the empty character vector for NaT
            %                   datetimes.
            %
            % See also YEAR, QUARTER, MONTH, WEEK, JULIANDATE, YMD.
            import matlab.internal.datetime.getDayNames
            
            if nargin == 1
                kind = 1;
            else
                kind = getChoice(kind,{'DayOfMonth' 'DoM' 'DayOfWeek' 'DoW' 'DayOfYear' 'DoY' 'Name' 'LongName' 'ShortName'},[1 1 2 2 3 3 4 4 5]);
            end
            ucal = getDateFieldsStructure;
            fieldIDs = [ucal.DAY_OF_MONTH ucal.DAY_OF_WEEK ucal.DAY_OF_YEAR ucal.DAY_OF_WEEK ucal.DAY_OF_WEEK];
            d = matlab.internal.datetime.getDateFields(this.data,fieldIDs(kind),this.tz);
            
            if kind > 3
                if kind == 4
                    names = getDayNames('long', getDatetimeSettings('locale'));
                else % kind == 5
                    names = getDayNames('short', getDatetimeSettings('locale'));
                end
                names{end+1} = ''; % return empty character vector for NaT
                d(isnan(d)) = length(names);
                d = reshape(names(d),size(this.data));
            end
        end
        
        function h = hour(this)
            %HOUR Hour numbers of datetimes.
            %   H = HOUR(T) returns the hour numbers of the datetimes in T. H is an
            %   array the same size as T containing integer values from 0 to 23.
            %
            % See also MINUTE, SECOND, TIMEOFDAY, HMS.
            ucal = getDateFieldsStructure;
            h = matlab.internal.datetime.getDateFields(this.data,ucal.HOUR_OF_DAY,this.tz);
        end
        
        function m = minute(this)
            %MINUTE Minute numbers of datetimes.
            %   M = MINUTE(T) returns the minute numbers of the datetimes in T. M is an
            %   array the same size as T containing integer values from 0 to 23.
            %
            % See also HOUR, SECOND, TIMEOFDAY, HMS.
            ucal = getDateFieldsStructure;
            m = matlab.internal.datetime.getDateFields(this.data,ucal.MINUTE,this.tz);
        end
        
        function s = second(this,kind)
            %SECOND Second numbers of datetimes.
            %   S = SECOND(T) returns the second values of the datetimes in T. S is an array
            %   the same size as T containing values (including a fractional part) from 0 to
            %   strictly less than 60.
            %
            %   For datetimes whose time zone is 'UTCLeapSeconds', SECOND returns a value
            %   from 60 to 61 for datetimes that are during a leap second occurrence.
            %
            %   S = SECOND(T,KIND) returns the kind of second values specified by KIND. KIND
            %   is one of the character vectors
            %
            %      secondofminute - The second of the minute. This is the default.
            %      secondofday    - The second of the day, from 0 to 86399.
            %
            % See also HOUR, MINUTE, TIMEOFDAY, HMS.
            if nargin == 1
                kind = 1;
            else
                kind = getChoice(kind,{'SecondOfMinute' 'SoM' 'SecondOfDay' 'SoD'},[1 1 2 2]);
            end
            ucal = getDateFieldsStructure;
            if kind == 1 % second+fraction within current minute
                s = matlab.internal.datetime.getDateFields(this.data,ucal.SECOND,this.tz);
            else         % second+fraction within current day
                s = matlab.internal.datetime.getDateFields(this.data,ucal.MILLISECOND_OF_DAY,this.tz) / 1000;
            end
        end
        
        function d = timeofday(this)
            %TIMEOFDAY Elapsed time since midnight for datetimes.
            %   D = TIMEOFDAY(T) returns an array of durations equal to the elapsed time since
            %   midnight for each of the datetimes in T, i.e. T - DATESHIFT(T,'START','DAY').
            %
            %   For unzoned datetimes, and in most other cases, D is equal to
            %
            %      HOURS(T.Hour) + MINUTES(T.Minute) + SECONDS(T.Second)
            %
            %   However, for datetimes whose TimeZone property is set to a time zone that observes
            %   Daylight Saving Time, on days where a Daylight Saving Time shift occurs, for times
            %   after the shift occurs, D differs from that sum by the amount of the shift.
            %
            %   Examples:
            %
            %      Create an unzoned datetime array, get the hour, minute, and second
            %      properties, and compare to the elapsed time since midnight:
            %         t = datetime(2015,3,8) + hours(1:4)
            %         [hrs,mins,secs] = hms(t)
            %         d = timeofday(t)
            %
            %      Set the times of day in one unzoned datetime array according to the
            %      times of day in another unzoned datetime array:
            %         t1 = datetime(2015,3,7) + hours(1:4)
            %         t2 = datetime(2015,3,repmat(8,1,4))
            %         t2 = dateshift(t2,'start','day') + timeofday(t1)
            %
            %      Create a zoned datetime array, on a day with a Daylight Saving Time
            %      shift, get the hour, minute, and second properties, and compare to the
            %      elapsed time since midnight:
            %         tz = 'America/New_York';
            %         fmt = 'dd-MMM-yyyy HH:mm:ss z';
            %         t = datetime(2015,3,8,'TimeZone',tz,'Format',fmt) + hours(1:4)
            %         [hrs,mins,secs] = hms(t)
            %         d = timeofday(t)
            %
            %      Set the times of day in one datetime array according to the times of day in
            %      another datetime array. This method works regardless of the time zone or the
            %      day of year. In 'America/New_York', 2:00AM did not exist on 8-Mar-2015, and
            %      so that element becomes 3:00AM:
            %         t1 = datetime(2015,3,7) + hours(1:4)
            %         tz = 'America/New_York';
            %         fmt = 'dd-MMM-yyyy HH:mm:ss z';
            %         t2 = datetime(2015,3,repmat(8,1,4),'TimeZone',tz,'Format',fmt)
            %         t2.Hour = t1.Hour; t2.Minute = t1.Minute; t2.Second = t1.Second
            %
            %   See also HMS, HOUR, MINUTE, SECOND, DATESHIFT, DURATION.
            d = this - dateshift(this,'start','day');
        end
        
        function [tz,dst] = tzoffset(this)
            %TZOFFSET Time zone offset of datetimes.
            %   DT = TZOFFSET(T) returns an array of durations equal to the time zone offset
            %   from UTC of the datetimes in T. For datetimes that occur in Daylight Saving
            %   Time, DT includes the time shift for DST. In other words, DT is the amount
            %   of time that each datetime in T differs from UTC.
            %
            %   The offset for unzoned datetimes is not defined.
            %
            %   [DT,DST] = TZOFFSET(T) also returns the time shift for Daylight Saving Time
            %   for each datetime in T.
            %
            %   See also TIMEZONE, ISDST.
            if isempty(this.tz)
                tz = NaN(size(this.data));
                dst = tz;
            else
                ucal = getDateFieldsStructure;
                [tz,dst] = matlab.internal.datetime.getDateFields(this.data,[ucal.ZONE_OFFSET ucal.DST_OFFSET],this.tz);
            end
            % Add the raw offset and the DST offset to get the total offset
            tz = duration(0,0,tz+dst,'Format','hh:mm');
            if nargout > 1
                dst = duration(0,0,dst,'Format','hh:mm');
            end
        end
        
        % no need for eomday method, that's day(dateshift(t,'end','month'))
        % no need for weekday method, that's day(t,'dayofweek')
        
        % These return logicals
        
        function tf = isweekend(this)
            %ISWEEKEND True for datetimes occurring on a weekend.
            %   TF = ISWEEKEND(T) returns a logical vector the same size as the datetime
            %   array T, with logical 1 (true) in elements where the corresponding element
            %   of T is a datetime the occurs on a weekend day, and logical 0 (false)
            %   otherwise.
            %
            %   See also ISDST, DAY.
            ucal = getDateFieldsStructure;
            dow = matlab.internal.datetime.getDateFields(this.data,ucal.DAY_OF_WEEK,this.tz);
            tf = (dow == 1) | (dow == 7); % Sunday is always day 1, regardless of locale
        end
        
        function tf = isdst(this)
            %ISDST True for datetimes occurring during Daylight Saving Time.
            %   TF = ISDST(T) returns a logical vector the same size as the datetime array
            %   T, with logical 1 (true) in elements where the corresponding element of T is
            %   a datetime the occurs during Daylight Saving Time, and logical 0 (false)
            %   otherwise.
            %
            %   ISDST returns false for datetimes in an "unzoned" array.
            %
            %   See also ISWEEKEND, TZOFFSET, TIMEZONE.
            ucal = getDateFieldsStructure;
            tf = matlab.internal.datetime.getDateFields(this.data,ucal.DST_OFFSET,this.tz) ~= 0;
            tf(isnan(this.data)) = false;
        end
        
        
        %% Array methods
        function [varargout] = size(this,varargin)
            [varargout{1:nargout}] = size(this.data,varargin{:});
        end
        function l = length(this)
            l = length(this.data);
        end
        function n = ndims(this)
            n = ndims(this.data);
        end
        
        function n = numel(this,varargin)
            if nargin == 1
                n = numel(this.data);
            else
                n = numel(this.data,varargin{:});
            end
        end
        
        function t = isempty(a),  t = isempty(a.data);  end
        function t = isscalar(a), t = isscalar(a.data); end
        function t = isvector(a), t = isvector(a.data); end
        function t = isrow(a),    t = isrow(a.data);    end
        function t = iscolumn(a), t = iscolumn(a.data); end
        function t = ismatrix(a), t = ismatrix(a.data); end
        
        function result = cat(dim,varargin)
            try
                [argsData,result] = datetime.isequalUtil(varargin);
            catch ME
                if strcmp(ME.identifier,'MATLAB:datetime:IncompatibleTZ')
                    error(message('MATLAB:datetime:cat:IncompatibleTZ'));
                elseif strcmp(ME.identifier,'MATLAB:datetime:InvalidComparison')
                    error(message('MATLAB:datetime:cat:InvalidConcatenation'));
                else
                    throw(ME);
                end
            end
            result.data = cat(dim,argsData{:}); % use fmt/tz from the first array
        end
        function result = horzcat(varargin)
            try
                result = cat(2,varargin{:});
            catch ME
                throw(ME);
            end
        end
        function result = vertcat(varargin)
            try
                result = cat(1,varargin{:});
            catch ME
                throw(ME);
            end
        end

        function this = ctranspose(this)
            try
                this.data = transpose(this.data); % NOT ctranspose
            catch ME
                throwAsCaller(ME);
            end
        end
        function this = transpose(this)
            try
                this.data = transpose(this.data);
            catch ME
                throwAsCaller(ME);
            end
        end
        function this = reshape(this,varargin)
            this.data = reshape(this.data,varargin{:});
        end
        function this = permute(this,order)
            this.data = permute(this.data,order);
        end
        
        %% Relational operators
        function t = eq(a,b)
            %EQ Equality comparison for datetimes.
            %   A == B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where the elements are equal, logical 0 (false)
            %   otherwise.  A and B must be the same size, or either can be a scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %  
            %   C = EQ(A,B) is called for the syntax 'A == B'.
            %
            %   See also NE, LT, LE, GE, GT.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) == 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = ne(a,b)
            %NE Not-equality comparison for datetimes.
            %   A ~= B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where the elements are not equal, and logical 0
            %   (false) otherwise.  A and B must be the same size, or either can be a
            %   scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %  
            %   C = NE(A,B) is called for the syntax 'A ~= B'.
            %
            %   See also EQ, LT, LE, GE, GT.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) ~= 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = lt(a,b)
            %LT Less than or equal comparison for datetimes.
            %   A < B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where A(I) < B(I), and logical 0 (false) otherwise.
            %   A and B must be the same size, or either can be a scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %
            %   C = LT(A,B) is called for the syntax 'A < B'.
            %
            %   See also EQ, NE, LE, GE, GT.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) < 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = le(a,b)
            %LE Less than or equal comparison for datetimes.
            %   A <= B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where A(I) <= B(I), and logical 0 (false) otherwise.
            %   A and B must be the same size, or either can be a scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %
            %   C = LE(A,B) is called for the syntax 'A <= B'.
            %
            %   See also EQ, NE, LT, GE, GT.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) <= 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = ge(a,b)
            %GE Greater than or equal comparison for datetimes.
            %   A >= B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where A(I) >= B(I), and logical 0 (false) otherwise.
            %   A and B must be the same size, or either can be a scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %
            %   C = GE(A,B) is called for the syntax 'A >= B'.
            %
            %   See also EQ, NE, LT, LE, GT.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) >= 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = gt(a,b)
            %GT Greater than comparison for datetimes.
            %   A > B returns a logical matrix the same size as the datetime arrays A and B
            %   with logical 1 (true) where A(I) > B(I), and logical 0 (false) otherwise.  A
            %   and B must be the same size, or either can be a scalar.
            %
            %   A or B can also be a date string, an array of date strings, a date character vector, 
            %   or a cell array of date character vectors.
            %
            %   C = GT(A,B) is called for the syntax 'A > B'.
            %
            %   See also EQ, NE, LT, LE, GE.
            try
                [aData,bData] = datetime.compareUtil(a,b);
                t = relopSign(aData,bData) > 0;
            catch ME
                throwAsCaller(ME);
            end
        end
        
        function t = isequal(varargin)
            %ISEQUAL True if datetime arrays are equal.
            %   TF = ISEQUAL(A,B) returns logical 1 (true) if the datetime arrays A and B
            %   are the same size and contain the same values, and logical 0 (false)
            %   otherwise. Either A or B may also be a cell array of character vectors, 
            %   or an array of strings.
            %
            %   TF = ISEQUAL(A,B,C,...) returns logical 1 (true) if all the input arguments
            %   are equal.
            %
            %   NaT elements are not considered equal to each other.  Use ISEQUALN to treat
            %   NaT elements as equal.
            %
            %   See also ISEQUALN, EQ.
            narginchk(2,Inf);
            try
                argsData = datetime.isequalUtil(varargin);
            catch ME
                if strcmp(ME.identifier,'MATLAB:datetime:InvalidComparison')
                    % silently return false
                elseif any(strcmp(ME.identifier,{'MATLAB:datetime:IncompatibleTZ' 'MATLAB:datetime:IncompatibleTZLeapSeconds'}))
                    % silently return false
                elseif any(strcmp(ME.identifier,{'MATLAB:datetime:AutoConvertString' 'MATLAB:datetime:AutoConvertStrings'}))
                    warning(message('MATLAB:datetime:AutoConvertStrings'));
                else
                    throw(ME);
                end
                t = false;
                return
            end
            t = isequal(argsData{:});
        end
        
        function t = isequaln(varargin)
            %ISEQUALN True if datetime arrays are equal, treating NaT elements as equal.
            %   TF = ISEQUALN(A,B) returns logical 1 (true) if the datetime arrays A and B
            %   are the same size and contain the same values or corresponding NaT elements,
            %   and logical 0 (false) otherwise.  Either A or B may also be a cell array of
            %   character vectors, or an array of strings.
            %
            %   TF = ISEQUALN(A,B,C,...) returns logical 1 (true) if all the input arguments
            %   are equal.
            %
            %   Use ISEQUAL to treat NaT elements as unequal.
            %
            %   See also ISEQUAL, EQ.
            narginchk(2,Inf);
            try
                argsData = datetime.isequalUtil(varargin);
            catch ME
                if strcmp(ME.identifier,'MATLAB:datetime:InvalidComparison')
                    % silently return false
                elseif any(strcmp(ME.identifier,{'MATLAB:datetime:IncompatibleTZ' 'MATLAB:datetime:IncompatibleTZLeapSeconds'}))
                    % silently return false
                elseif any(strcmp(ME.identifier,{'MATLAB:datetime:AutoConvertString' 'MATLAB:datetime:AutoConvertStrings'}))
                    warning(message('MATLAB:datetime:AutoConvertStrings'));
                else
                    throw(ME);
                end
                t = false;
                return
            end
            t = isequaln(argsData{:});
        end
        
    end % public methods block
    
    methods(Hidden = true)
        %% Display and strings
        function disp(this,name)
            if (nargin < 2)
                name = '';
            end
            
            try
                displayFun(this,name);
            catch ME
                throwAsCaller(ME);
            end
        end

        %% Arrayness
        function n = end(this,k,n)
            try
                n = builtin('end',this.data,k,n);
            catch ME
                throwAsCaller(ME);
            end
        end
        
        %% Subscripting
        this = subsasgn(this,s,rhs)
        that = subsref(this,s)

        function sz = numArgumentsFromSubscript(~,~,~)
            % This function is for internal use only and will change in a
            % future release.  Do not use this function.
            sz = 1;
        end

        %% Variable Editor methods
        % These functions are for internal use only and will change in a
        % future release.  Do not use this function.
        [out,warnmsg] = variableEditorClearDataCode(this, varname, rowIntervals, colIntervals)
        [out,warnmsg] = variableEditorColumnDeleteCode(this, varName, colIntervals)
        out = variableEditorInsert(this, orientation, row, col, data)
        metadata = variableEditorMetadata(this)
        [out,warnmsg] = variableEditorMetadataCode(this, varName, index, propertyName, propertyString)
        out = variableEditorPaste(this, rows, columns, data)
        [out,warnmsg] = variableEditorRowDeleteCode(this, varName, rowIntervals)
        [out,warnmsg] = variableEditorSetDataCode(this, varname, row, col, rhs)
        [out,warnmsg] = variableEditorSortCode(~, varName, columnIndexStrings, direction)
        
        %% Error stubs
        % Methods to override functions and throw helpful errors
        function d = double(d), error(message('MATLAB:datetime:InvalidNumericConversion','double')); end %#ok<MANU>
        function d = single(d), error(message('MATLAB:datetime:InvalidNumericConversion','single')); end %#ok<MANU>
        function d = months(d), error(message('MATLAB:datetime:NoMonthsMethod')); end %#ok<MANU>
        function d = timezone(d), error(message('MATLAB:datetime:NoTimeZoneMethod')); end %#ok<MANU>
        function d = format(d), error(message('MATLAB:datetime:NoFormatMethod')); end %#ok<MANU>
        function d = floor(varargin), error(message('MATLAB:datetime:UseDateshiftMethod','floor')); end %#ok<STOUT>
        function d = ceil(varargin), error(message('MATLAB:datetime:UseDateshiftMethod','ceil')); end %#ok<STOUT>
        function d = round(varargin), error(message('MATLAB:datetime:UseDateshiftMethod','round')); end %#ok<STOUT>
    end % hidden public methods block
    
    methods(Hidden = true, Static = true)
        function d = empty(varargin)
        %EMPTY Create an empty datetime array.
        %   D = DATETIME.EMPTY() creates a 0x0 datetime array.
        %
        %   D = DATETIME.EMPTY(M,N,...) or D = DATETIME.EMPTY([N M ...]) creates
        %   an N-by-M-by-... datetime array.  At least one of N,M,... must be zero.
        %
        %   See also DATETIME.
            if nargin == 0
                d = datetime([],[],[]);
            else
                dData = zeros(varargin{:});
                if numel(dData) ~= 0
                    error(message('MATLAB:datetime:empty:InvalidSize'));
                end
                d = datetime([],[],[]);
                d.data = dData;
            end
        end
        
        function tzs = allTimeZones()
        % This function is for internal use only and will change in a
        % future release.  Do not use this function.
        tzs = cell2table(matlab.internal.datetime.getDefaults('TimeZones'), ...
            'VariableNames',{'Name' 'CanonicalName' 'UTCOffset' 'DSTOffset'});
        end
        
        function t = createInternal(data,format)
        % This function is for internal use only and will change in a
        % future release.  Do not use this function.
        %
        %   T = CREATEINTERNAL(DATA, FORMAT)
        %   
        %   Create a datetime with the internal datatype DATA and input
        %   format FORMAT.
        %
            t = datetime;
            t.data = data;
            t.tz = '';
            t.fmt = format; 
        end
        
        function [allfmts,numUnambiguous] = getFormatsForGuessing(locale)
        % This function is for internal use only and will change in a
        % future release.  Do not use this function.

            % This function rarely needs to change. To add formats to this
            % function, one must be congnizant of all the callers of this
            % function which may/may not need to know about the new formats
            % being introduced.
            %
            % One has to also be wary of the M/d d/M ambiguity and
            % carefully select where the formats to be added must be
            % placed.
        
            % Generate a list of DATE formats which are common and completely
            % distinguishable from each other. These are also paired with TIME
            % formats to generate a list of datetime formats to use when
            % guessing.
            datefmts = {'dd-MMM-uuuu' ... % MATLAB default
                        'MMM/dd/uuuu' ... % Leading month name. Does not accept numbers for month.
                        'uuuu-MM-dd' ...  % ISO8601, has month number but is sufficiently recognizable
                        matlab.internal.datetime.getDefaults('LocaleFormat',locale,'uuuuMMMd') % locale's MMM default
                       };

            % Timestamp formats which are completely distinguishable from each
            % other. HH can take values from 0:23, hh can take values from
            % 1:12 and requires AM/PM to be correctly recognized. Also check
            % for '' which means no time stamp.
            timefmts = {''  ' HH:mm:ss' ' HH:mm' ' hh:mm:ss aa',' hh:mm aa'};


            dtfmts = strcat(repmat(datefmts',1,length(timefmts)),repmat(timefmts,length(datefmts),1));
            % These formats are not paired with times. e.g. 'MMM-uuuu HH:mm' will not
            % be guessed.
            dtfmts = [dtfmts(:); 'MMM-uuuu';'QQQ-uuuu';'HH:mm:ss';'hh:mm:ss aa';'uuuuMMdd''T''HHmmss'];

            % It's possible that the locale-format matches another format in
            % the list.
            dtfmts = unique(dtfmts,'stable');
            
            % Get the preferred MMdd format for the locale.  This may start
            % with MM/dd or dd/MM. If that is the case, we want to take both
            % cases MM/dd... and dd/MM... 
            mmddfmts = {matlab.internal.datetime.getDefaults('LocaleFormat',locale,'uuuuMMdd')};
            if ~isempty(regexp(mmddfmts{1},'^MM/dd','once'))
                % Create a day-first equivalent for this month-first locale
                mmddfmts{2} = regexprep(mmddfmts{1},'MM/dd','dd/MM');
            elseif ~isempty(regexp(mmddfmts{1},'^dd/MM','once'))
                % Create a month-first equivalent for this day-first locale
                mmddfmts{2} = regexprep(mmddfmts{1},'dd/MM','MM/dd');
            end
            % If there was only one MMdd format used for this locale, return
            % the total length as the number of unambiguous, otherwise, just
            % use the number of dtfmts.
            numUnambiguous = numel(dtfmts) + numel(timefmts)*isscalar(mmddfmts);
            % Generate all the pairs of date and time formats.
            mdfmts = strcat(repmat(mmddfmts',1,length(timefmts)),repmat(timefmts,length(mmddfmts),1));
            
            allfmts = [dtfmts(:); mdfmts(:)];
        end
        
        function fmt = getDefaultFormatForLocale(locale)
        % This function is for internal use only and will change in a
        % future release.  Do not use this function.
            if nargin == 0
                locale = matlab.internal.datetime.getDefaults('locale');
            else
                locale = matlab.internal.datetime.verifyLocale(locale);
            end
            fmt = matlab.internal.datetime.getDefaults('localeformat',locale,'uuuuMMdd');
        end
        
        function setLocalTimeZone(tz)
        % This function is for internal use only and will change in a
        % future release.  Do not use this function.
            if nargin == 0, tz = []; end
            sessionLocalTimeZone(tz);
        end
    end % hidden static public methods block
    
    methods(Static, Access='public')
        setDefaultFormats(format,formatStr)
    end
    
    %% Protected methods
    methods(Access='protected')
        displayFun(this,objectname)
        this = subsasgnDot(this,s,rhs)
        this = subsasgnParens(this,s,rhs)
        value = subsrefDot(this,s)
        value = subsrefParens(this,s)
        fmt = getDisplayFormat(this)
    end % protected methods block

    methods(Static, Access='protected')
        [a,b,prototype] = compareUtil(a,b)
        [a,b] = arithUtil(a,b)
        [args,prototype] = isequalUtil(args)
        t = convertFrom(value,type,tz,epoch)
    end % static protected methods block
end % classdef


%%%%%%%%%%%%%%%%% Local functions %%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------
function tf = isDateOnlyFormat(fmt)
tf = isempty(regexp(fmt,'[hHkKmsSaA]','once'));
end

%-----------------------------------------------------------------------
function pivot = verifyPivot(pivot)

import matlab.internal.tableUtils.isScalarInt

if ~isScalarInt(pivot)
    error(message('MATLAB:datetime:InvalidPivotYear'));
end
end

%-----------------------------------------------------------------------
function thisData = currentTime(tz)
try %#ok<ALIGN>
    % Get time components from the system clock and create an unzoned internal
    % value to match those components. clock returns seconds as d.p. rounded to
    % the nearest milli- (Windows) or micro-second (Linux/OSX). datetime counts
    % in ms, so (re)rounding 1000*sec to 3 digits (microsec) gets the first case
    % exactly right and more or less leaves the second case alone. For the second
    % case, could round in datetime's higher precision, but the cosmetic tweak
    % later on would break that anyway, so don't bother.
    c = clock;
    thisData = matlab.internal.datetime.createFromDateVec({c(1) c(2) c(3) c(4) c(5) 0 round(1000*c(6),3)},'');
    
    % The assumption is that the system clock is not leap-second-aware. In the
    % unlikely event that it _is_, and returns a time during a 60th second,
    % the above shifts it into the 0th second of the next minute. But if the
    % requested result is in 'UTCLeapSeconds', will need to preserve the leap
    % second below.
    
    if nargin == 0 || isempty(tz)
        % If there's a session value set to override the system time zone, shift the
        % system clock's result as if the system was set to the specified time zone.
        sessionLocalTZ = sessionLocalTimeZone();
        if ~isempty(sessionLocalTZ)
            thisData = timeZoneAdjustment(thisData,sessionLocalTZ,'');
            thisData = timeZoneAdjustment(thisData,'',datetime.SystemTimeZone);
        end
        tz = '';
    else
        % The requested result is zoned, remove the system time zone's offset (raw
        % offset plus DST) from the system clock to translate the value to UTC. This
        % leaves the actual time unchanged, but in general (unless the requested tz
        % is the same as the system tz), it gives a different timestamp.
        thisData = timeZoneAdjustment(thisData,'',datetime.SystemTimeZone);
        
        if strcmp(tz,datetime.UTCLeapSecsZoneID)
            % The internal value was created without accounting for leap seconds, add them.
            thisData = addLeapSeconds(thisData);
            if c(6) >= 60 % system returned a current time during a leap second
                % Shift the time from the 0th second back into the 60th second.
                thisData = matlab.internal.datetime.datetimeSubtract(thisData,1000,true);
            end
        end
    end
    
    % For cosmetics, make (today + (now - today)) equal to now.
    ucal = getDateFieldsStructure;
    t0 = matlab.internal.datetime.datetimeFloor(thisData,ucal.DAY_OF_MONTH,tz);
    dt = matlab.internal.datetime.datetimeSubtract(thisData,t0);
    thisData = matlab.internal.datetime.datetimeAdd(t0,dt);
catch ME, throwAsCaller(ME), end
end

%-----------------------------------------------------------------------
function inData = expandNumericInputs(inData)
sz = [1 1];
for i = 1:length(inData)
    field = inData{i};
    if ~isscalar(field)
        sz = size(field);
        break
    end
end
for i = 1:length(inData)
    field = inData{i};
    if isscalar(field)
        inData{i} = repmat(full(double(field)),sz);
    else % let createFromDateVec check for size mismatch
        inData{i} = full(double(field));
    end
end
end

%-----------------------------------------------------------------------
function handleParseErrors(inData,supplied,fmt,inputFmt,locale)

try %#ok<ALIGN>
    if isempty(locale)
        locale = matlab.internal.datetime.getDefaults('locale');
    end
    if supplied.InputFormat || supplied.Format
        if supplied.InputFormat, fmt = inputFmt; end
        if supplied.Locale
            if isscalar(inData)
                error(message('MATLAB:datetime:ParseErrWithLocale',inData{1},fmt,locale));
            else
                error(message('MATLAB:datetime:ParseErrsWithLocale',fmt,locale));
            end
        elseif ~isempty(regexp(fmt,'[eMz]{3,}','match','once'))
            if isscalar(inData)
                error(message('MATLAB:datetime:ParseErrSuggestLocale',inData{1},fmt,locale));
            else
                error(message('MATLAB:datetime:ParseErrsSuggestLocale',fmt,locale));
            end
        else
            if isscalar(inData)
                error(message('MATLAB:datetime:ParseErr',inData{1},fmt));
            else
                error(message('MATLAB:datetime:ParseErrs',fmt));
            end
        end
    else % guessing a format
        if supplied.Locale
            if isscalar(inData)
                error(message('MATLAB:datetime:UnrecognizedDateStringWithLocale',inData{1},locale));
            else
                error(message('MATLAB:datetime:UnrecognizedDateStringsWithLocale',locale));
            end
        else
            if isscalar(inData)
                error(message('MATLAB:datetime:UnrecognizedDateStringSuggestLocale',inData{1},locale));
            else
                error(message('MATLAB:datetime:UnrecognizedDateStringsSuggestLocale',locale));
            end
        end
    end
    
catch ME, throwAsCaller(ME); end
end
