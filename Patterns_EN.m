function Stats = Patterns_EN(samples, ch, n, k, c)
% ----------------------------------------------
% Patterns_EN - Counting Patterns with Exactly k-2 Failures
%
% Usage:
%   Stats = Patterns_EN(samples, ch, n, k, c)
%
% Description:
%   Patterns_EN counts the number of occurrences of a pattern with exactly
%   k-2 failures in a matrix of sequences where each row represents a sample.
%
% Input:
%   - samples: Matrix where each row represents a sample.
%   - ch: Character to search for in the sequences.
%   - n: Length of the sequence.
%   - k: Number of required successes separated by exactly k-2 failures.
%   - c: Flag for circular sequences (0 for non-circular, 1 for circular).
%
% Output:
%   - Stats: Matrix containing the counts of occurrences of the specified
%            pattern for each sample.
%
% Example:
%   samples = [
%      1, 0, 1, 1, 2, 2, 1, 0, 1, 1, 2, 0, 0, 0, 1, 1, 1, 2, 1, 1;
%      0, 1, 1, 2, 1, 1, 0, 1, 1, 2, 0, 0, 1, 0, 1, 1, 1, 2, 0, 0
%   ];
%   character_to_search = 1;
%   sequence_length = 20;
%   required_successes = 3;
%   circular_flag = 0;
%   result = Patterns_EN(samples, character_to_search, sequence_length, required_successes, circular_flag);
%
% ----------------------------------------------

if nargin < 4
    error('Patterns:Patterns_EN:TooFewInputs', 'Input arguments are undefined.');
elseif nargin < 5
    c = 0;
end

[nrows, ncols] = size(samples);

if ncols < n
    error('Patterns:Patterns_EN:WrongInputs', 'Column size and sequence length are not matching');
elseif nrows == 1 && mod(ncols, n) == 0
    samples = vec2mat(samples, ncols);
end

r = length(k);
Stats = repmat(0, nrows, r);
for d = 1:r
    if c == 1
        tpsamples = [samples samples(:, [1:k(d) - 1])];
    else
        tpsamples = samples;
    end
    [nrows, ncols] = size(tpsamples);

    if k(d) > n
        no = 0;
    else
        for i = 1:nrows
            tline1 = tpsamples(i, :);
            sad = find(tline1 == ch);
            j = 1;
            no(i) = 0;

            while j <= length(sad) - 1
                if (sad(j + 1) - sad(j)) == (k(d) - 1)
                    no(i) = no(i) + 1;
                    j = j + 2;
                else
                    j = j + 1;
                end
            end
        end
    end
    clear temp j
    Stats(:, d) = no';
end
clear no tpsamples;
end

