function Stats = Patterns_AMN(samples, ch, n, k, c)
% PATTERNS_AMN Extracts patterns from a matrix of samples.
%
% STATS = PATTERNS_AMN(SAMPLES, CH, N, K, C) analyzes the input matrix
% SAMPLES and counts the occurrences of patterns with a specified
% character CH and a maximum separation of K-2 failures. The function returns
% the statistics in the STATS matrix.
%
% Input Arguments:
%   - SAMPLES: Matrix of samples where patterns are extracted.
%   - CH: Character to be considered in pattern extraction.
%   - N: Sequence length.
%   - K: Maximum separation allowed between two successes.
%   - C: Circular extension flag (optional, default is 0).
%
% Output Argument:
%   - STATS: Matrix containing the count of patterns for each specified K.

% Check the number of input arguments
if nargin < 4
    error('Patterns:Patterns_AMN:TooFewInputs', 'Input arguments are undefined.');
elseif nargin < 5
    c = 0;
end

% Check the size of the input matrix
[nrows, ncols] = size(samples);
if ncols < n
    error('Patterns:Patterns_AMN:WrongInputs', 'Column size and sequence length do not match.');
elseif nrows == 1 && mod(ncols, n) == 0
    samples = vec2mat(samples, ncols);
end

% Initialize the Stats matrix
r = length(k);
Stats = repmat(0, nrows, r);

% Loop through each pattern length in K
for d = 1:r
    % Extend samples if circular flag is set
    if c == 1
        tpsamples = [samples samples(:, 1:k(d) - 1)];
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
                if (sad(j + 1) - sad(j)) <= (k(d) - 1)
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
