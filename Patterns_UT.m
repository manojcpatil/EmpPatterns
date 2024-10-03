function Stats = Patterns_UT(samples, userpattern, r, c)
% ----------------------------------------------
% ----------------------------------------------
% Patterns_UT - Counting Non-overlapping Patterns
%
% Usage:
%   Stats = Patterns_UT(samples, userpattern, r, c)
%
% Description:
%   Patterns_UT counts the number of trials required to observe a specified number
%   (r) of non-overlapping occurrences of a user-defined pattern in a matrix of 
%   sequences where each row represents a sample.
%
% Input:
%   - samples: Matrix where each row represents a sample.
%   - userpattern: Pattern to search for in the sequences. Should be a vector.
%   - r: Number of non-overlapping occurrences of the specified pattern to count.
%   - c: Flag for circular sequences (0 for non-circular, 1 for circular).
%
% Output:
%   - Stats: Matrix containing the number of trials required to observe r
%            non-overlapping occurrences of the specified patterns for each sample.
%
% Example:
%   samples = [
%      1, 0, 1, 1, 2, 2, 1, 0, 1, 1, 2, 0, 0, 0, 1, 1, 1, 2, 1, 1;
%      0, 1, 1, 2, 1, 1, 0, 1, 1, 2, 0, 0, 1, 0, 1, 1, 1, 2, 0, 0
%   ];
%   pattern = [1, 2];  % Pattern is now a vector
%   r = 2;  % Number of non-overlapping occurrences to count
%   circular_flag = 0;
%   result = Patterns_UT(samples, pattern, r, circular_flag);
%
%
% ----------------------------------------------

if nargin<2
    error('Patterns:Patterns_UN:TooFewInputs','Input arguments are undefined.');
elseif nargin<3
    c=0;
end

[nrows,n]=size(samples);
Stats = repmat(n+1,nrows, 1);
k=length(userpattern);

if nrows==1
    samples=vec2mat(samples,n);
end

for d=1:length(k)
    if c==1
        tpsamples=[samples samples(:,[1:k(d)-1])];
    else
        tpsamples=samples;
    end
    [nrows,ncols]=size(tpsamples);

    if k(d)>n
        no=0;
    else

        for i=1:nrows
            ncols_temp=ncols;
            temp = tpsamples(i,:);
            j=k(d);
            no(i)=0;
            rth=0;

            while j<=ncols_temp && rth<r
                if min(temp(j-k(d)+1:j)==userpattern)==1
                    rth=rth+1;
                    if rth==r
                        no(i)=j;
                        break;
                    end
                    j=j+k(d);                    
                else
                    j=j+1;
                end
            end            
        end
        clear temp j ncols_temp
    end
    Stats(:,d)=no';

    clear no tpsamples;
end
end