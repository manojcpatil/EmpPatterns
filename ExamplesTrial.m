samples=[3,1,2,2,3,2,2,2,2,0,3,2,3,1,3,3,0,2,3,2
2,3,2,2,3,1,2,2,3,2,2,3,2,1,3,1,2,1,2,1
0,0,3,2,2,3,2,1,3,2,1,2,2,3,2,0,2,2,2,3
2,1,2,3,2,1,2,0,1,2,1,3,1,2,2,2,2,3,2,1
3,0,1,1,1,2,3,2,1,2,2,2,2,3,1,2,3,2,1,3]

% Parameters
ch = 2; % Character to be considered
n = 20; % Sequence length
k = 3;  % Pattern length
r = 2;  % Occurrence number
c = 0;  % Circular extension flag (optional, default is 0)

N1 = Patterns_AMN(samples, ch, n, k, c);
T1 = Patterns_AMT(samples, ch, k, r, c);
N2 = Patterns_EN(samples, ch, n, k, c);
T2 = Patterns_ET(samples, ch, k, r, c);
N3 = Patterns_ALN(samples, ch, n, k, c);
T3 = Patterns_ALT(samples, ch, k, r, c);

userpattern = [2, 3, 2] % Example user-defined pattern
Nu = Patterns_UN(samples, userpattern, c);
Tu = Patterns_UT(samples, userpattern, r, c);

disp('[N1,T1,N2,T2,N3,T3,Nu,Tu]')
[N1,T1,N2,T2,N3,T3,Nu,Tu]


% Parameters
ch = 2; % Character to be considered
n = 20; % Sequence length
k = 3;  % Pattern length
r = 2;  % Occurrence number
c = 0;  % Circular extension flag (optional, default is 0)

M1 = Patterns_AMM(samples, ch, n, k, c);
W1 = Patterns_AMW(samples, ch, k, r, c);
M2 = Patterns_EM(samples, ch, n, k, c);
W2 = Patterns_EW(samples, ch, k, r, c);
M3 = Patterns_ALM(samples, ch, n, k, c);
W3 = Patterns_ALW(samples, ch, k, r, c);

userpattern = [2, 3, 2] % Example user-defined pattern
Mu = Patterns_UM(samples, userpattern, c);
Wu = Patterns_UW(samples, userpattern, r, c);

[M1,W1,M2,W2,M3,W3,Mu,Wu]

samples=binornd(3,0.6,5,20)
ch = 2; n = 20; k = [3,4,5]; c = 0; r = 2;
JN1 = Patterns_AMN(samples, ch, n, k, c)
JT1 = Patterns_AMT(samples, ch, k, r)
[JN1,JT1]