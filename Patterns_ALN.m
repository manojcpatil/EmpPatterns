function Stats = Patterns_ALN(samples, ch, n, k, c)
% PATTERNS_ALN Extracts patterns from a matrix of samples.
%
% STATS = PATTERNS_ALN(SAMPLES, CH, N, k, C) analyzes the input matrix
% SAMPLES and counts the occurrences of patterns with a specified
% character CH and various pattern lengths specified by the array PATTERNLENGTHS.
% The function returns the statistics in the STATS matrix.
%
% Input Arguments:
%   - SAMPLES: Matrix of samples where patterns are extracted.
%   - CH: Character to be considered in pattern extraction.
%   - N: Sequence length.
%   - k: Array of pattern lengths to analyze.
%   - C: Circular extension flag (optional, default is 0, non-circular).
%
% Output Argument:
%   - STATS: Matrix containing the count of patterns for each pattern length.
%
% Example:
%   samples = [1 2 1 1 2 1; 2 1 2 1 2 1];
%   ch = 1;
%   n = 6;
%   k = [2, 3];
%   c = 0;
%   stats = Patterns_ALN(samples, ch, n, k, c);
%
% See also: VEC2MAT
%
% Author: Manoj C Patil
% Date: Today's Date
%
% Version History:
%   1.0 - Initial version.
%


% Check the number of input arguments
if nargin < 4
    error('Patterns:Patterns_ALN:TooFewInputs', 'Input arguments are undefined.');
elseif nargin < 5
    c = 0;
end

% Check the size of the input matrix
[nrows, ncols] = size(samples);
if ncols < n
    error('Patterns:Patterns_ALN:WrongInputs', 'Column size and sequence length do not match.');
elseif nrows == 1 && mod(ncols, n) == 0
    samples = vec2mat(samples, ncols);
end

% Initialize the Stats matrix
r = length(k);
Stats = repmat(0, nrows, r);

% Loop through each pattern length in k
for d = 1:r
    % Extend samples if circular flag is set
    if c == 1
        tpsamples = [samples samples(:, [1:k(d) - 1])];
    else
        tpsamples = samples;
    end

    % Get the size of the extended samples matrix
    [nrows, ncols] = size(tpsamples);

    % Count the occurrences of patterns for each row
    if k(d) > n
        no = 0;
    else
        for i = 1:nrows
            tline1 = tpsamples(i, :);
            sad = find(tline1 == ch);
            j = 1;
            no(i) = 0;

            while j <= length(sad) - 1
                if (sad(j + 1) - sad(j)) > k(d) - 2
                    no(i) = no(i) + 1;
                    j = j + 2;
                else
                    j = j + 1;
                end
            end
        end
    end

    % Store the count in the Stats matrix
    Stats(:, d) = no';
end

% Clear temporary variables
clear no tpsamples;

end
