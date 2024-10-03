function Stats = Patterns_AMM(samples, ch, n, k, c)
% PATTERNS_AMM Extracts patterns from a matrix of samples.
%
% STATS = PATTERNS_AMM(SAMPLES, CH, N, K, C) analyzes the input matrix
% SAMPLES and counts the occurrences of patterns with a specified
% character CH and a minimum length K. The function returns the
% statistics in the STATS matrix.
%
% Input Arguments:
%   - SAMPLES: Matrix of samples where patterns are extracted.
%   - CH: Character to be considered in pattern extraction.
%   - N: Sequence length.
%   - K: Array of pattern lengths to analyze.
%   - C: Circular extension flag (optional, default is 0).
%
% Output Argument:
%   - STATS: Matrix containing the count of patterns for each length K.
%
% Example:
%   samples = [1 2 1 1 2 1; 2 1 2 1 2 1];
%   ch = 1;
%   n = 6;
%   k = [2, 3];
%   c = 0;
%   stats = Patterns_AMM(samples, ch, n, k, c);
%
% Author: Your Name
% Date: Today's Date
%
% Version History:
%   1.0 - Initial version.
%
% Note: Add any additional information, examples, or explanations as needed.

% Input validation
if nargin < 4
    error('Patterns:Patterns_AMM:TooFewInputs', 'Input arguments are undefined.');
elseif nargin < 5
    c = 0;
end

% Check the size of the input matrix
[nrows, ncols] = size(samples);
if ncols < n
    error('Patterns:Patterns_AMM:WrongInputs', 'Column size and sequence length do not match.');
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
    [nrows, n] = size(tpsamples);

    % Count the occurrences of patterns for each row
    if k(d) <= n
        no = zeros(nrows, 1);
        for i = 1:nrows
            tline1 = tpsamples(i, :);
            sad = find(tline1 == ch);
            j = 1;
            no(i) = 0;

            while j <= length(sad) - 1
                if (sad(j + 1) - sad(j)) <= (k(d) - 1)
                    no(i) = no(i) + 1;
                end
                j = j + 1;
            end
        end

        % Store the count in the Stats matrix
        Stats(:, d) = no;
    end
end

end
