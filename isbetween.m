function boolean = isbetween(a, b, c, varargin)
% Is "a" between values "b" and "c"?
%
% boolean = isbetween(a, b, c, ...
%                    'inclusive', inclusive_switch, ...
%                    'threshold', threshold_value])
%
% Example: isbetween(1, 0, 3) returns true.
%          isbetween(1:4, 0, 3) returns [true, true, true, false]
%          isbetween(1:4, 0, 3, false) returns [true, true, false, false]
%
% Inputs:
%   a, b, c             scalar or vector values
%
% Options:
%   inclusive_switch        Boolean: if a equals b or c, is it "between"?
%                               default: true (equals counts as "between")
%                               otherwise, specify false
%   threshold_value         Is there some small value for which overlaps
%                               should be ignored?
% Output:
%   boolean             true if a is between b & c

% Kevin J. Delaney
% November 19, 2008
% January 25, 2010      Extended to "b", "c" being vectors.
% March 10, 2010        Allow "a" to be a vector.
% April 5, 2010         Added 'inclusive', 'threshold' options.
% May 20, 2010          Return an output same orientation as input.

boolean = [];
inclusive = true;
threshold_value = 0;

if ~exist('a', 'var')
    help(mfilename);
    return
end

if ~isnumeric(a)
	errordlg('Input "a" is non-numeric.', mfilename);
    return
end

if ~exist('b', 'var') || ~isnumeric(b)
	errordlg('Input "b" is missing or non-numeric.', mfilename);
    return
end

if ~exist('c', 'var') || ~isnumeric(c)
	errordlg('Input "c" is missing or non-numeric.', mfilename);
    return
end

for index = 1 : 2 : (length(varargin) - 1)
    option_name = varargin{index};
    
    if isempty(option_name) || ~ischar(option_name)
        errordlg('Option name is empty or non-char.', mfilename);
        return
    end
    
    option_value = varargin{index + 1};
    
    if isempty(option_value)
        errordlg('Option value is empty.', mfilename);
        return
    end
    
    switch lower(option_name)
                   
        case 'inclusive'

            if ischar(option_value)
                switch lower(option_value)
                    case {'inclusive', 'yes', 'true', 'y', 't'}
                        inclusive = true;
                    case {'exclusive', 'false', 'no', 'f', 'n'}
                        inclusive = false;
                    otherwise
                        errordlg(['Unknown "inclusive" selection "', option_value, '".'], ...
                            mfilename);
                        return
                end
            elseif islogical(option_value)
                inclusive = option_value;
            else
                errordlg('Option "inclusive" neither boolean nor char.', mfilename);
                return
            end

                        
        case 'threshold'
            
            if ~isnumeric(option_value)
                errordlg('Input "threshold" is non-numeric.', mfilename);
                return
            end
            
            threshold_value = option_value;
            
            
        otherwise
            errordlg(['Unknown option "', option_name, '".'], mfilename);
            return
    end
end

%   Are b & c vectors?
num_tests = length(b);

if length(c) ~= num_tests
    errordlg(['Length mismatch. There are ', num2str(length(b)), ...
        ' "b" values but ', num2str(length(c)), ' "c" values.'], mfilename);
    return
end

%   Make sure both are column vectors.
b = reshape(b, [num_tests, 1]);
c = reshape(c, [num_tests, 1]);

lower_limit = min([b, c], [], 2);
upper_limit = max([b, c], [], 2);

% if num_dims(a) > 1
%     errordlg('Sorry--can''t accomodate "a" with more than one dimension.', ...
%         mfilename);
%     return
% end

boolean = false(num_tests, length(a));

for a_index = 1:length(a)
    if inclusive
        boolean(:, a_index) = (a(a_index) >= (lower_limit + threshold_value)) & (a(a_index) <= (upper_limit - threshold_value));
    else
        boolean(:, a_index) = (a(a_index) >  (lower_limit + threshold_value)) & (a(a_index) <  (upper_limit - threshold_value));
    end
end

%   Return an output same orientation as input.
if size(a, 1) > size(a, 2)
    boolean = boolean';
end