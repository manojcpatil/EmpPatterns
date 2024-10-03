function Stats = Patterns_EM(samples, ch, n, k, c)
if nargin<4
    error('Patterns:Patterns_ALM:TooFewInputs','Input arguments are undefined.');
elseif nargin<5
    c=0;
end

[nrows,ncols]=size(samples);

if ncols<n
    error('Patterns:Patterns_ALM:WrongInputs','Column size and sequence length are not matching');
elseif nrows==1 && mod(ncols,n)==0
    samples=vec2mat(samples,ncols);
end

r=length(k);
Stats=repmat(0,nrows,r);
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
            tline1 = tpsamples(i,:);
            sad = find(tline1==ch);
            j=1;
            no(i)=0;

            while j <= length(sad) - 1
                if (sad(j + 1) - sad(j)) == (k(d) - 1)
                    no(i) = no(i) + 1;
                end
                j = j + 1;
            end
        end
    end
    clear temp j
    Stats(:,d)=no';
end
clear no tpsamples;
end