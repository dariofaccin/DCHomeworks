clc; close all; clear global; clearvars;

L = 63;             % length of PN sequence
Nh = 6;             % Bound on the length of h

sigdB = -8;         % noise variance
sigmaw = 10^(sigdB/10);
w = wgn(4*L,1,sigdB);
w_0 = w([1:2:end-1]);          % take even samples for h0
w_1 = w([2:2:end]);            % take odd samples for h1

PN = PNSeq(L);      % ML sequence 
x=[PN ; PN];

a1 = -0.9635;       % IIR filter given
a2 = 0.4642;
h = impz(1, [1 a1 a2]);
[h_even,h_odd] = polyphase(h,Nh);    % Polyphase decomposition

%% ESTIMATE OF h0 WITH THE CORRELATION METHOD
z_0 = filter(h_even, 1, x);       
d_0 = z_0 + w_0;
h0_cor = corr_method(x, d_0);       % correlation method

%% ESTIMATE OF h1 WITH THE CORRELATION METHOD
z_1=filter(h_odd, 1, x);
d_1 = z_1 + w_1;
h1_cor=corr_method(x, d_1);         % correlation method

if Nh<L
    h0_cor=h0_cor(1:ceil(Nh/2));
    h1_cor=h1_cor(1:floor(Nh/2));
end
h_cor = PS(h0_cor, h1_cor);


%% ESTIMATE WITH THE LS METHOD
h0_ls=LS(x, d_0, L);     % least-square methods
h1_ls=LS(x, d_1, L);

if Nh<L
    h0_ls=h0_ls(1:ceil(Nh/2));
    h1_ls=h1_ls(1:floor(Nh/2));
end
h_ls = PS(h0_ls, h1_ls);
h_ls=h_ls(1:Nh);

%% VALUES FOR THE TABLE
table = zeros(Nh,3)
for i=1:Nh
    table(i,:) = [h(i) h_cor(i) h_ls(i)];
end
table
















