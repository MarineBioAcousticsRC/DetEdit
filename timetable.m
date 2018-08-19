classdef (Sealed) timetable < tabular
%TIMETABLE Timetable.
%   Timetables are used to collect heterogeneous data and metadata into a single
%   container, where each row has a timestamp.  Timetables are suitable for
%   storing column-oriented or tabular timestamped data that are often stored as
%   columns in a text file or in a spreadsheet.  Timetables can accommodate
%   variables of different types, sizes, units, etc.  They are often used to
%   store experimental data, with rows representing different observations and
%   columns representing different measured variables.
%
%   Use the TIMETABLE constructor to create a timetable from variables in the
%   MATLAB workspace.  Use TABLE2TIMETABLE or ARRAY2TIMETABLE to create a
%   timetable from a table or array, respectively. Use the readtable function to
%   create a table by reading data from a text or spreadsheet file, and then use
%   TABLE2TIMETABLE to convert that to a timetable.
%
%   Timetables can be subscripted using parentheses much like ordinary numeric
%   arrays, but in addition to numeric and logical indices, you can use a
%   timetable's variable and row names as indices.  You can access individual
%   variables in a timetable much like fields in a structure, using dot
%   subscripting.  You can access the contents of one or more variables using
%   brace subscripting.
%
%   Timetables can contain different kinds of variables, including numeric,
%   logical, character, categorical, and cell.  However, a timetable is a
%   different class than the variables that it contains.  For example, even a
%   timetable that contains only variables that are double arrays cannot be
%   operated on as if it were itself a double array.  However, using dot
%   subscripting, you can operate on a variable in a timetable as if it were a
%   workspace variable.  Using brace subscripting, you can operate on one or
%   more variables in a timetable as if they were in a homogeneous array.
%
%   A timetable TT has properties that store metadata such as its variable names
%   and row times.  Access or assign to a property using P = TT.Properties.PropName
%   or TT.Properties.PropName = P, where PropName is one of the following:
%
%   TIMETABLE metadata properties:
%       Description          - A character vector describing the timetable
%       DimensionNames       - A two-element cell array of character vectors containing names of
%                              the dimensions of the timetable
%       VariableNames        - A cell array containing names of the variables in the timetable
%       VariableDescriptions - A cell array of character vectors containing descriptions of the
%                              variables in the timetable
%       VariableUnits        - A cell array of character vectors containing units for the variables
%                              in the timetable
%       RowTimes             - A datetime or duration vector containing the times associated 
%                              with each row in the timetable. The time are not required to be
%                              uniformly-spaced, unique, or sorted, and RowTimes may contain
%                              missing values.
%       UserData             - A variable containing any additional information associated
%                              with the timetable.  You can assign any value to this property.
%
%   TIMETABLE methods and functions:
%     Construction and conversion:
%       timetable        - Create a timetable from workspace variables.
%       array2timetable  - Convert homogeneous array to timetable.
%       table2timetable  - Convert table to timetable.
%       timetable2table  - Convert timetable to table.
%     Size and shape:
%       istimetable      - True for timetables.
%       size             - Size of a timetable.
%       width            - Number of variables in a timetable.
%       height           - Number of rows in a timetable.
%       ndims            - Number of dimensions of a timetable.
%       numel            - Number of elements in a timetable.
%       horzcat          - Horizontal concatenation for timetables.
%       vertcat          - Vertical concatenation for timetables.
%     Set membership:
%       intersect        - Find rows common to two timetables.
%       ismember         - Find rows in one timetable that occur in another timetable.
%       setdiff          - Find rows that occur in one timetable but not in another.
%       setxor           - Find rows that occur in one or the other of two timetables, but not both.
%       unique           - Find unique rows in a timetable.
%       union            - Find rows that occur in either of two timetables.
%       join             - Merge two timetables by matching up rows using key variables.
%       innerjoin        - Inner join between two timetables.
%       outerjoin        - Outer join between two timetables.
%     Data manipulation and reorganization
%       retime           - Adjust a timetable and its data to a new time vector.
%       synchronize      - Synchronize timetables.
%       sortrows         - Sort rows of a timetable.
%       lag              - Lag or lead data in a timetable.
%       issorted         - TRUE for a sorted timeable.
%       isregular        - TRUE for a regular timetable.
%       summary          - Print summary of a timetable.
%       sortrows         - Sort rows of a timetable.
%       stack            - Stack up data from multiple variables into a single variable.
%       unstack          - Unstack data from a single variable into multiple variables.
%       ismissing        - Find elements in a timetable that contain missing values.
%       standardizeMissing - Insert missing data indicators into a timetable.
%     Computations on timetables:
%       varfun           - Apply a function to variables in a timetable.
%       rowfun           - Apply a function to rows of a timetable.
%     Subscripting into timetables:
%       timerange        - Timetable row subscripting by time range.
%       withtol          - Timetable row subscripting by time with tolerance.
%       vartype          - Timetable variable subscripting by variable type.
%
%   See also TIMETABLE, TABLE2TIMETABLE, ARRAY2TIMETABLE

%   Copyright 2016 The MathWorks, Inc.
        
    properties(Constant, GetAccess='protected')
        propertyNames = [fieldnames(tabular.arrayPropsDflts); ...
                                    matlab.internal.table.tableMetaDim.propertyNames; ...
                                    matlab.internal.table.tableVarNamesDim.propertyNames; ...
                                    matlab.internal.table.tableRowTimesDim.propertyNames];
        defaultDimNames = dfltTimetableDimNames();
    end
        
    properties(Transient, Access='protected')
        data = cell(1,0);
        
        % Create the metaDim to display the row labels header
        metaDim = matlab.internal.table.tableMetaDim(2,timetable.defaultDimNames,true);
        rowDim = matlab.internal.table.tableRowTimesDim(0,datetime.empty(0,1));
        varDim = matlab.internal.table.tableVarNamesDim(0);
        
        % 'Properties' will appear to contain this, as well as the per-row, per-var,
        % and per-dimension properties contained in rowDim, varDim. and metaDim,
        arrayProps = timetable.arrayPropsDflts;
    end
            
    methods
        function t = timetable(varargin)
        %TIMETABLE Create a timetable from workspace variables.
        %   Use TIMETABLE to create a timetable from variables in the MATLAB workspace.
        %   Use TABLE2TIMETABLE or ARRAY2TIMETABLE to create a timetable from a table or
        %   array, repectively. To create a timetable from data from a text or
        %   spreadsheet file, use READTABLE and then use TABLE2TIMETABLE to convert the
        %   result from a table to a timetable.
        %
        %   TT = TIMETABLE(ROWTIMES, VAR1, VAR2, ...) creates a timetable TT from the
        %   workspace variables VAR1, VAR2, ..., using the datetime or duration vector
        %   ROWTIMES as the time vector.  All variables must have the same number of
        %   rows.
        %
        %   TT = TIMETABLE(VAR1, VAR2, ..., 'RowTimes',ROWTIMES) creates a timetable
        %   using the specified datetime or duration vector, ROWTIMES, as the time
        %   vector. Other datetime or duration inputs become variables in TT.
        %
        %   TT = TIMETABLE(..., 'VariableNames', {'name1', ..., 'name_M'}) creates a
        %   timetable using the specified variable names. The names must be valid MATLAB
        %   identifiers, and unique.
        %
        %   Timetables can contain variables that are built-in types, or objects that
        %   are arrays and support standard MATLAB parenthesis indexing of the form
        %   var(i,...), where i is a numeric or logical vector that corresponds to rows
        %   of the variable. In addition, the array must implement a SIZE method with a
        %   DIM argument, and a VERTCAT method.
        %
        % See also TABLE2TIMETABLE, ARRAY2TIMETABLE.
        
            if nargin == 0
                % Nothing to do
            else
                % Construct from separate variables
                [numVars,numRows] = tabular.countVarInputs(varargin);
                vars = varargin(1:numVars);

                if numVars < nargin
                    pnames = {'VariableNames'  'RowTimes'};
                    dflts =  {            []           []};
                    try
                        [varnames,rowtimes,supplied] ...
                            = matlab.internal.table.parseArgs(pnames, dflts, varargin{numVars+1:end});
                    catch ME
                        % The inputs included a 1xM string that was interpreted as a param name,
                        % but something went wrong. If all of the preceedng inputs had one row,
                        % or if that string was first, these two errors suggest that the string
                        % was intended as data. Suggest alternative options, in that case.
                        if (numVars == 0 || numRows == 1)
                            if strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName') || ...
                                    strcmp(ME.identifier,'MATLAB:table:parseArgs:WrongNumberArgs')
                                cause = message('MATLAB:table:ConstructingFromStrings');
                                ME = addCause(ME,MException(cause.Identifier,'%s',getString(cause)));
                            end
                        end
                        throw(ME);
                    end
                else
                    supplied.VariableNames = false;
                    supplied.RowTimes = false;
                end

                if ~supplied.VariableNames
                    % Get the workspace names of the input arguments from inputname if
                    % variable names were not provided. Need these names before looking
                    % through vars for the time vector.
                    varnames = repmat({''},1,numVars);
                    for i = 1:numVars, varnames{i} = inputname(i); end
                end

                if supplied.RowTimes
                    if numVars == 0, numRows = length(rowtimes); end
                else
                    rowtimes = vars{1};
                    if ~isdatetime(rowtimes) && ~isduration(rowtimes)
                        if numVars == 1 && matlab.internal.datatypes.istabular(vars{1})
                            error(message('MATLAB:timetable:NoTimeVectorTableInput'));
                        else
                            error(message('MATLAB:timetable:NoTimeVector'));
                        end
                    end
                    vars(1) = [];
                    numVars = numVars - 1; % don't count time index in vars

                    if ~supplied.VariableNames
                        rowtimesName = varnames{1};
                        varnames(1) = [];
                        if ~isempty(rowtimesName)
                            % If the rows times came from a data input (not from the
                            % RowTimes param), and var names were not provided, get the
                            % row dim name from the inputs. Otherwise, leave the default
                            % row dim name alone.
                            t.metaDim = t.metaDim.setLabels(rowtimesName,1);
                        end
                    end
                end

                if ~supplied.VariableNames
                    % Fill in default names for data args where inputname couldn't. Do
                    % this after removing the time vector from the other vars, to get the
                    % default names numbered correctly.
                    empties = cellfun('isempty',varnames);
                    if any(empties)
                        varnames(empties) = matlab.internal.table.dfltVarNames(find(empties)); %#ok<FNDSB>
                    end
                    % Make sure default names or names from inputname don't conflict
                    varnames = matlab.lang.makeUniqueStrings(varnames,{},namelengthmax);
                end

                t.rowDim = matlab.internal.table.tableRowTimesDim(numRows,rowtimes);
                t.varDim = matlab.internal.table.tableVarNamesDim(numVars,varnames); % error if invalid, duplicate, or empty
                t.data = vars;
                
                % Detect conflicts between the var names and the default dim names.
                t.metaDim = t.metaDim.checkAgainstVarLabels(varnames);
            end
        end
    end
    
    methods(Hidden, Static)
        function t = empty(varargin)
        %EMPTY Create an empty table.
        %   T = TIMETABLE.EMPTY() creates a 0x0 timetable.
        %
        %   T = TIMETABLE.EMPTY(NROWS,NVARS) or T = TIMETABLE.EMPTY([NROWS NVARS]) creates
        %   an NROWSxNVARS timetable.  At least one of NROWS or NVARS must be zero. If
        %   NROWS is positive, T's time vector contains NaTs.
        %
        %   See also TIMETABLE, ISEMPTY.
            if nargin == 0
                t = timetable();
            else
                sizeOut = size(zeros(varargin{:}));
                if prod(sizeOut) ~= 0
                    error(message('MATLAB:timetable:empty:EmptyMustBeZero'));
                elseif length(sizeOut) > 2
                    error(message('MATLAB:timetable:empty:EmptyMustBeTwoDims'));
                else
                    % Create a 0x0 timetable, and then resize to the correct number
                    % of rows or variables.
                    t = timetable();
                    if sizeOut(1) > 0
                        t.rowDim = t.rowDim.lengthenTo(sizeOut(1));
                    end
                    if sizeOut(2) > 0
                        t.varDim = t.varDim.lengthenTo(sizeOut(2));
                        t.data = cell(1,sizeOut(2)); % assume double
                    end
                end
            end
        end
        
        % Called by cell2timetable, struct2timetable
        function t = fromScalarStruct(s)
            % This function is for internal use only and will change in a
            % future release.  Do not use this function.
            
            % Construct a timetable from a scalar struct
            vnames = fieldnames(s);
            p = length(vnames);
            if p > 0
                n = unique(structfun(@(f)size(f,1),s));
                if ~isscalar(n)
                    error(message('MATLAB:table:UnequalFieldLengths'));
                end
            else
                n = 0;
            end
            t = timetable();
            t.data = struct2cell(s)';
            t.rowDim = matlab.internal.table.tableRowTimesDim(n);
            t.varDim = matlab.internal.table.tableVarNamesDim(p,vnames);
        end
    end % hidden static methods block
    
    methods(Access = 'protected')  
        % Helper methods used in join/innerjoin/outerjoin
        [a_out,b_out,supplied,keys,leftKeys,rightKeys,leftVars,rightVars,mergeKeys] ...
            = prejoin(a,b,supplied,type,keys,leftKeys,rightKeys,leftVars,rightVars,mergeKeys);
        c = postjoin(c,a,b,type,leftVarDim)
        
        % Used by join functions to prepare timetable for tabular join
        t = copyRowTimesAsVariable(tt);
        
        function b = cloneAsEmpty(a)            
        %CLONEASEMPTY Create a new empty table from an existing one.
%             if strcmp(class(a),'timetable') %#ok<STISA>
                b = timetable();
%             else % b is a subclass of timetable
%                 b = a; % respect the subclass
%                 % leave b.metaDim alone;
%                 b.rowDim = matlab.internal.table.tableRowTimesDim();
%                 b.varDim = matlab.internal.table.tableVarNamesDim();
%                 b.data = cell(1,0);
%                 leave b.arrayProps alone
%             end
        end
        
        function errID = throwSubclassSpecificError(obj,msgid,varargin)
            % Throw the timetable version of the msgid error, using varargin as the
            % variables to fill the holes in the message.
            errID = throwSubclassSpecificError@tabular(obj,['timetable:' msgid],varargin{:});
            if nargout == 0
                throwAsCaller(errID);
            end
        end
    end % protected methods block

    %%%% PERSISTENCE BLOCK ensures correct save/load across releases %%%%%%
    %%%% Properties and methods in this block maintain the exact class %%%%
    %%%% schema required for TIMETABLE to persist through MATLAB releases %    
    properties(Constant, Access='private')
        % running timetable version. This is used only for managing forward
        % compatibility. Value is not saved when an instance is serialized
        version = 2.0;
    end
    
    methods(Hidden)
        function tt_serialized = saveobj(tt)
            tt_serialized = struct;
            tt_serialized.versionSavedFrom = tt.version;  % scalar double. version number of saved timetable
            tt_serialized.minCompatibleVersion = 2.0;     % scalar double. minimum running version required to reconstruct a timetable from this serialized data
            tt_serialized.arrayProps = tt.arrayProps;     % scalar struct. Fields hold array specific properties not derived from the above
            tt_serialized.data       = tt.data;           % 1-by-numVars cell vector. Each cell corresponds to data from one variable
            tt_serialized.numDims    = tt.metaDim.length; % scalar double. Number of dimensions
            tt_serialized.dimNames   = tt.metaDim.labels; % 1-by-numDims cell of char vector. Names of each dimension
            tt_serialized.numRows    = tt.rowDim.length;  % scalar double. Number of rows
            tt_serialized.rowTimes   = tt.rowDim.labels;  % numRows-by-1 datetime or duration (must exist). Time of each row of data
            tt_serialized.numVars    = tt.varDim.length;  % scalar double. Number of variables
            tt_serialized.varNames   = tt.varDim.labels;  % 1-by-numVars cell of char vector. Names, if any, of each variable
        end
    end
    
    methods(Hidden, Static)
        function tt = loadobj(tt_serialized)
            % Always default construct an empty timetable, and recreate a
            % proper timetable in the current schema using attributes
            % loaded from the serialized struct            
            tt = timetable();
            
            % Return an empty instance if current version is below the
            % minimum compatible version of the serialized object
            if timetable.version < tt_serialized.minCompatibleVersion
                warning(message('MATLAB:timetable:IncompatibleLoad'));
                return;
            end
            
            % Restore core data and dimension objects
            if isstruct(tt_serialized.rowTimes) % Optimized rowTimes is saved, expand into a full time vector in this version
                % An optimized rowTimes struct should have these fields:
                %     origin  : scalar datetime or duration. Start of the regularly spaced time vector
                %     stepSize: scalar duration. Time step of the regularly spaced time vector
                rowTimes = tt_serialized.rowTimes.origin + ( 0:tt_serialized.rowTimes.stepSize:tt_serialized.rowTimes.stepSize*(tt_serialized.numRows-1) );
            else
                rowTimes = tt_serialized.rowTimes;
            end
            tt.rowDim  = matlab.internal.table.tableRowTimesDim(tt_serialized.numRows, rowTimes);
            tt.varDim  = matlab.internal.table.tableVarNamesDim(tt_serialized.numVars, tt_serialized.varNames);
            tt.metaDim = matlab.internal.table.tableMetaDim    (tt_serialized.numDims, tt_serialized.dimNames, true);
            tt.data    = tt_serialized.data;
            
            % Restore arrayProps struct from compatible data in the
            % serialized object. If data was saved from a different
            % version, remove any fields in the serialized arrayProps
            % that are not recognized in the current version.
            if tt_serialized.versionSavedFrom ~= timetable.version
                toRemove = setdiff(fieldnames(tt_serialized.arrayProps),fieldnames(tabular.arrayPropsDflts));
                if ~isempty(toRemove)
                    tt_serialized.arrayProps = rmfield(tt_serialized.arrayProps,toRemove);
                end
            end
            tt = tt.setProperties(tt_serialized.arrayProps); % Restore the remaining metadata properties.
        end
    end
end

%-------------------------------------------------------------------------------
function names = dfltTimetableDimNames()
names = { getString(message('MATLAB:timetable:uistrings:DfltRowDimName')) ...
          getString(message('MATLAB:timetable:uistrings:DfltVarDimName')) };
end
