clc; close all; clear global; clearvars;

L=127;
Nlim=20;

% Additive noise
sigdB=-8;
sigmaw=10^(sigdB/10);
w_0=sigmaw*randn(2*L,1);
w_1=sigmaw*randn(2*L,1);

%% ML sequence
% The input is the repetition of two ML sequences, each of length L
x=[PNSeq(L); PNSeq(L)];

%% POLYPHASE REALIZATION
% First we compute the impulse response, then we split in even and
% odd samples (polyphase components)
a1 = -0.9635;
a2 = 0.4642;
h = impz(1, [1 a1 a2]);
h = h(1:Nlim);

% Even samples
h_even=zeros(Nlim/2,1);
for k=1:(Nlim/2)
    h_even(k)=h(2*k-1);
end

% Odd samples
h_odd=zeros(Nlim/2,1);
for k=1:Nlim/2
    h_odd(k)=h(2*k);
end

%% ESTIMATE OF h WITH THE CORRELATION METHOD

r_0=filter(h_even, 1, x);
r_1=filter(h_odd, 1, x);

h0_est=r_dx(x, r_0);
h1_est=r_dx(x, r_1);

% Estimated impulse response: take even/odd samples from the polyphase
% components (S/P)

h_est=zeros(Nlim,1);
for i=1:Nlim
    if (i<=L)
    h_est(2*i-1)=h0_est(i);
    h_est(2*i)=h1_est(i);
    end
end

h_est=h_est(1:Nlim);
% Plot analytical and estimated (via correlation method) impulse response
figure('Name','Analytical and estimated (Correlation method) impulse response');
stem(0:19,h_est);
hold on;
stem(0:19,h,'r*');
title('h_{analytic} vs h_{estimate-CORR}');
xlabel('n');
ylim([-0.5 1.2]); xlim([-2 20]);
legend('h_{est-CORR}','h_{analytic}')

%% LS
d0=r_0;
d1=r_1;

h0_ls=LS(x, d0, L);
h1_ls=LS(x, d1, L);

h_ls=zeros(Nlim,1);
for i=1:Nlim
    if (i<=L)
    h_ls(2*i-1)=h0_ls(i);
    h_ls(2*i)=h1_ls(i);
    end
end
h_ls=h_ls(1:Nlim);
figure('Name','Analytical and estimated (LS method) impulse response');
stem(0:19, h_ls);
hold on;
stem(0:19,h,'r*');
title('h_{analytic} vs h_{estimate-LS}');
xlabel('n');
ylim([-0.5 1.2]); xlim([-2 20]);
legend('h_{est-LS}','h_{analytic}');

%% sigma_w ESTIMATION