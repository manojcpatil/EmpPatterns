function Stats = Patterns_AMW(samples, ch, k, r, c)
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
            tpsamples = [samples samples(:, [1:k(d)-1])];
        else
            tpsamples = samples;
        end

        % Initialize a vector to store results for each row
        no = repmat(ncols+1,nrows, 1);

        % Loop through each row of samples
        for i = 1:nrows
            tline1 = tpsamples(i, :);
            sad = find(tline1 == ch);

            % Search for patterns in the row
            temp = 0;
            j = 1;
            while j <= (length(sad) - 1) && temp < r
                if (sad(j + 1) - sad(j)) <= (k(d) - 1)
                    temp = temp + 1;                
                    if temp == r
                        no(i) = sad(j + 1);
                        break;
                    end
                end
                j = j + 1;
            end
        end

        % Store results in the Stats cell array
        Stats(:, d) = no;
    end
end
