function Stats = Patterns_ALT(samples, ch, k, r, c)
% PATTERNS_ALT Extracts patterns from a matrix of samples.
%
% STATS = PATTERNS_ALT(SAMPLES, CH, K, R, C) analyzes the input matrix
% SAMPLES and finds the positions of the Rth occurrence of a pattern
% with a specified character CH and a minimum length K. The function
% returns the positions in the STATS matrix.
%
% Input Arguments:
%   - SAMPLES: Matrix of samples where patterns are extracted.
%   - CH: Character to be considered in pattern extraction.
%   - K: Array of pattern lengths to analyze.
%   - R: The occurrence number to find for each pattern length.
%   - C: Circular shift flag (optional, default is 0).
%
% Output Argument:
%   - STATS: Matrix containing the positions of the Rth occurrences
%            for each pattern length K. If no such occurrence is found,
%            the corresponding entry is set to ncols + 1.
%
% Example:
%   samples = [1 2 1 1 2 1; 2 1 2 1 2 1];
%   ch = 1;
%   k = [2, 3];
%   r = 2;
%   c = 0;
%   stats = Patterns_ALT(samples, ch, k, r, c);
%
% Author: Manoj C. Patil
% Date: Today's Date
%
% Version History:
%   1.0 - Initial version.
%
% Note: Add any additional information, examples, or explanations as needed.

% Input validation
if nargin < 4
    error('Patterns:Patterns_ALT:TooFewInputs', 'Input arguments are undefined.');
elseif nargin < 5
    c = 0;
end

% Convert single-row input to a matrix
[nrows, ncols] = size(samples);
if nrows == 1
    samples = vec2mat(samples, ncols);
end

% Initialize Stats matrix as a cell array
    Stats = repmat(ncols+1,nrows, length(k));

% Loop through each pattern length
for d = 1:length(k)
    % Circular shift samples if c is 1
    if c == 1
        tpsamples = [samples samples(:, [1:k(d) - 1])];
    else
        tpsamples = samples;
    end

    % Initialize a vector to store results for each row
    no = repmat(ncols + 1, nrows, 1);

    % Loop through each row of samples
    for i = 1:nrows
        tline1 = tpsamples(i, :);
        sad = find(tline1 == ch);

        % Search for patterns in the row
        temp = 0;
        j = 1;
        while j <= (length(sad) - 1) && temp < r
            if (sad(j + 1) - sad(j)) > (k(d) - 2)
                temp = temp + 1;
                if temp == r
                    no(i) = sad(j + 1);
                    break;
                end
                j = j + 2;
            else
                j = j + 1;
            end
        end
    end

    % Store results in the Stats cell array
    Stats(:, d) = no;
end

end
