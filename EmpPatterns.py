import numpy as np

function Stats = Patterns_ALM(samples, ch, n, k, c)
    if nargin < 4
        error('Patterns:Patterns_ALM:TooFewInputs', 'Input arguments are undefined.');
    elseif nargin < 5
        c = 0;
    end

    [nrows, ncols] = size(samples);

    if ncols < n
        error('Patterns:Patterns_ALM:WrongInputs', 'Column size and sequence length do not match.');
    elseif nrows == 1 && mod(ncols, n) == 0
        samples = vec2mat(samples, ncols);
    end

    r = length(k);
    Stats = zeros(nrows, r);

    for d = 1:r
        if c == 1
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape
        no = np.zeros(nrows)

        if k[d] <= n:
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                no[i] = 0

                while j <= len(sad) - 2:
                    if (sad[j + 1] - sad[j]) == (k[d] - 1):
                        no[i] += 1
                    j += 1

            Stats[:, d] = no

    return Stats



def Patterns_ALN(samples, ch, n, k, c=0):
    # Check the size of the input matrix
    nrows, ncols = samples.shape

    if ncols < n:
        raise ValueError('Column size and sequence length do not match.')

    if nrows == 1 and ncols % n == 0:
        samples = np.reshape(samples, (ncols // n, n))

    # Initialize the Stats matrix
    r = len(k)
    Stats = np.zeros((nrows, r))

    # Loop through each pattern length in k
    for d in range(r):
        # Extend samples if circular flag is set
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Get the size of the extended samples matrix
        nrows, ncols = tpsamples.shape

        # Count the occurrences of patterns for each row
        if k[d] > n:
            no = np.zeros(nrows)
        else:
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                no = np.zeros(nrows)

                while j <= len(sad) - 2:
                    if (sad[j + 1] - sad[j]) > k[d] - 2:
                        no[i] += 1
                        j += 2
                    else:
                        j += 1

        # Store the count in the Stats matrix
        Stats[:, d] = no

    return Stats





def Patterns_ALT(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a numpy array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) > (k[d] - 2):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                    j += 2
                else:
                    j += 1

        # Store results in the Stats numpy array
        Stats[:, d] = no

    return Stats

def Patterns_ALW(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a numpy array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) > (k[d] - 2):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                j += 1

        # Store results in the Stats numpy array
        Stats[:, d] = no

    return Stats





def Patterns_AMM(samples, ch, n, k, c=0):
    # Input validation
    if c is None:
        c = 0

    # Check the size of the input matrix
    nrows, ncols = samples.shape
    if ncols < n:
        raise ValueError('Column size and sequence length do not match.')

    if nrows == 1 and ncols % n == 0:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize the Stats matrix
    r = len(k)
    Stats = np.zeros((nrows, r))

    # Loop through each pattern length in k
    for d in range(r):
        # Extend samples if circular flag is set
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Get the size of the extended samples matrix
        nrows, ncols = tpsamples.shape

        # Count the occurrences of patterns for each row
        if k[d] <= ncols:
            no = np.zeros(nrows)
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                while j <= len(sad) - 2:
                    if (sad[j + 1] - sad[j]) <= (k[d] - 1):
                        no[i] += 1
                    j += 1

            # Store the count in the Stats matrix
            Stats[:, d] = no

    return Stats

def Patterns_AMN(samples, ch, n, k, c=0):
    # Input validation
    if c is None:
        c = 0

    # Check the size of the input matrix
    nrows, ncols = samples.shape
    if ncols < n:
        raise ValueError('Column size and sequence length do not match.')

    if nrows == 1 and ncols % n == 0:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize the Stats matrix
    r = len(k)
    Stats = np.zeros((nrows, r))

    # Loop through each pattern length in k
    for d in range(r):
        # Extend samples if circular flag is set
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Get the size of the extended samples matrix
        nrows, ncols = tpsamples.shape

        # Count the occurrences of patterns for each row
        if k[d] > n:
            no = np.zeros(nrows)
        else:
            no = np.zeros(nrows)
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                while j <= len(sad) - 2:
                    if (sad[j + 1] - sad[j]) <= (k[d] - 1):
                        no[i] += 1
                        j += 2
                    else:
                        j += 1

            # Store the count in the Stats matrix
            Stats[:, d] = no

    return Stats






def Patterns_AMT(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a cell array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) <= (k[d] - 1):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                    j += 2
                else:
                    j += 1

        # Store results in the Stats cell array
        Stats[:, d] = no

    return Stats

def Patterns_AMW(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a cell array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) <= (k[d] - 1):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                j += 1

        # Store results in the Stats cell array
        Stats[:, d] = no

    return Stats

def Patterns_EM(samples, ch, n, k, c=0):
    # Input validation
    if c is None:
        c = 0

    nrows, ncols = samples.shape

    if ncols < n:
        raise ValueError('Column size and sequence length are not matching')

    if nrows == 1 and ncols % n == 0:
        samples = np.reshape(samples, (ncols, 1))

    r = len(k)
    Stats = np.zeros((nrows, r))

    for d in range(r):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                no[i] = 0

                while j <= len(sad) - 1:
                    if (sad[j + 1] - sad[j]) == (k[d] - 1):
                        no[i] += 1
                    j += 1

        Stats[:, d] = no

    return Stats




def Patterns_EN(samples, ch, n, k, c=0):
    # Input validation
    if c is None:
        c = 0

    # Check the size of the input matrix
    nrows, ncols = samples.shape
    if ncols < n:
        raise ValueError('Column size and sequence length are not matching')

    if nrows == 1 and ncols % n == 0:
        samples = np.reshape(samples, (ncols, 1))

    r = len(k)
    Stats = np.zeros((nrows, r))

    for d in range(r):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)
            for i in range(nrows):
                tline1 = tpsamples[i, :]
                sad = np.where(tline1 == ch)[0]
                j = 0
                no[i] = 0

                while j <= len(sad) - 2:
                    if (sad[j + 1] - sad[j]) == (k[d] - 1):
                        no[i] += 1
                        j = j + 2
                    else:
                        j = j + 1

        Stats[:, d] = no

    return Stats

def Patterns_ET(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a cell array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) == (k[d] - 1):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                    j += 2
                else:
                    j += 1

        # Store results in the Stats cell array
        Stats[:, d] = no

    return Stats

def Patterns_EW(samples, ch, k, r, c=0):
    # Input validation
    if c is None:
        c = 0

    # Convert single-row input to a matrix
    nrows, ncols = samples.shape
    if nrows == 1:
        samples = np.reshape(samples, (ncols, 1))

    # Initialize Stats matrix as a cell array
    Stats = np.full((nrows, len(k)), ncols + 1)

    # Loop through each pattern length
    for d in range(len(k)):
        # Circular shift samples if c is 1
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        # Initialize a vector to store results for each row
        no = np.full(nrows, ncols + 1)

        # Loop through each row of samples
        for i in range(nrows):
            tline1 = tpsamples[i, :]
            sad = np.where(tline1 == ch)[0]

            # Search for patterns in the row
            temp = 0
            j = 0
            while j <= len(sad) - 2 and temp < r:
                if (sad[j + 1] - sad[j]) == (k[d] - 1):
                    temp += 1
                    if temp == r:
                        no[i] = sad[j + 1]
                        break
                j += 1

        # Store results in the Stats cell array
        Stats[:, d] = no

    return Stats

def Patterns_UM(samples, userpattern, c=0):
    # Input validation
    if c is None:
        c = 0

    nrows, n = samples.shape
    k = len(userpattern)

    if nrows == 1:
        samples = np.reshape(samples, (n, 1))

    r = len(k)
    Stats = np.zeros((nrows, r))

    for d in range(r):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)

            for i in range(nrows):
                ncols_temp = ncols
                temp = tpsamples[i, :]
                j = k[d]
                no[i] = 0

                while j <= ncols_temp:
                    if np.min(temp[j - k[d] + 1:j] == userpattern) == 1:
                        no[i] = no[i] + 1
                    j = j + 1

            Stats[:, d] = no

    return Stats



def Patterns_UN(samples, userpattern, c=0):
    # Input validation
    if c is None:
        c = 0

    nrows, n = samples.shape
    k = len(userpattern)

    if nrows == 1:
        samples = np.reshape(samples, (n, 1))

    r = len(k)
    Stats = np.zeros((nrows, r))

    for d in range(r):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)

            for i in range(nrows):
                ncols_temp = ncols
                temp = tpsamples[i, :]
                j = k[d]
                no[i] = 0

                while j <= ncols_temp:
                    if np.min(temp[j - k[d] + 1:j] == userpattern) == 1:
                        no[i] = no[i] + 1
                        j = j + k[d]
                    else:
                        j = j + 1

            Stats[:, d] = no

    return Stats

def Patterns_UT(samples, userpattern, r, c=0):
    # Input validation
    if c is None:
        c = 0

    nrows, n = samples.shape
    Stats = np.full((nrows, 1), n + 1)
    k = len(userpattern)

    if nrows == 1:
        samples = np.reshape(samples, (n, 1))

    for d in range(len(k)):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)

            for i in range(nrows):
                ncols_temp = ncols
                temp = tpsamples[i, :]
                j = k[d]
                no[i] = 0
                rth = 0

                while j <= ncols_temp and rth < r:
                    if np.min(temp[j - k[d] + 1:j] == userpattern) == 1:
                        rth = rth + 1
                        if rth == r:
                            no[i] = j
                            break
                        j = j + k[d]
                    else:
                        j = j + 1

            Stats[:, d] = no

    return Stats

def Patterns_UW(samples, userpattern, r, c=0):
    # Input validation
    if c is None:
        c = 0

    nrows, n = samples.shape
    k = len(userpattern)

    if nrows == 1:
        samples = np.reshape(samples, (n, 1))

    Stats = np.full((nrows, len(k)), n + 1)

    for d in range(len(k)):
        if c == 1:
            tpsamples = np.concatenate((samples, samples[:, :k[d] - 1]), axis=1)
        else:
            tpsamples = samples

        nrows, ncols = tpsamples.shape

        if k[d] > n:
            no = 0
        else:
            no = np.zeros(nrows)

            for i in range(nrows):
                ncols_temp = ncols
                temp = tpsamples[i, :]
                j = k[d]
                no[i] = 0
                rth = 0

                while j <= ncols_temp and rth < r:
                    if np.min(temp[j - k[d] + 1:j] == userpattern) == 1:
                        rth = rth + 1
                        if rth == r:
                            no[i] = j
                            break
                    j = j + 1

            Stats[:, d] = no

    return Stats

