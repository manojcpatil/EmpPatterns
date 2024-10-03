function Stats = Patterns_UN(samples, userpattern, c)
% ----------------------------------------------
% Patterns_UN - Counting Non-overlapping Patterns
%
% Usage:
%   Stats = Patterns_UN(samples, userpattern, c)
%
% Description:
%   Patterns_UN counts the number of non-overlapping occurrences of a user-defined
%   pattern in a matrix of sequences where each row represents a sample.
%
% Input:
%   - samples: Matrix where each row represents a sample.
%   - userpattern: Pattern to search for in the sequences. Should be a vector.
%   - c: Flag for circular sequences (0 for non-circular, 1 for circular).
%
% Output:
%   - Stats: Matrix containing the counts of non-overlapping occurrences of the
%            specified patterns for each sample.
%
% Example:
%   samples = [
%      1, 0, 1, 1, 2, 2, 1, 0, 1, 1, 2, 0, 0, 0, 1, 1, 1, 2, 1, 1;
%      0, 1, 1, 2, 1, 1, 0, 1, 1, 2, 0, 0, 1, 0, 1, 1, 1, 2, 0, 0
%   ];
%   pattern = [1, 2];  % Pattern is now a vector
%   circular_flag = 0;
%   result = Patterns_UN(samples, pattern, circular_flag);
%
%
% ----------------------------------------------


if nargin<2
    error('Patterns:Patterns_UN:TooFewInputs','Input arguments are undefined.');
elseif nargin<3
    c=0;
end

[nrows,n]=size(samples);
k=length(userpattern);

if nrows==1
    samples=vec2mat(samples,n);
end

r=length(k);
for d=1:r
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

            while j<=ncols_temp
                if min(temp(j-k(d)+1:j)==userpattern)==1
                    no(i)=no(i)+1;j=j+k(d);
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