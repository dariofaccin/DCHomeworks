clc; close all; clear global; clearvars;

% Input signal
load('rec_input.mat');
T = 1;
L = 1023;
sigma_a = 2;

snr_db = 7;
snr_lin = 10^(snr_db/10);

% Matched filter
gm = conj(qc(end:-1:1));

% h
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);

% sample h with period T
h_T = downsample(h,4);

% figure()
% subplot(121), stem(-8:8,h,'b');
% xlabel('nT/4'), xlim([-9 9]), grid on
% subplot(122), stem(-2:2,h_T,'Color','red');
% xlabel('nT'), xlim([-3 3]), grid on

r_c_prime = filter(gm,1,r_c);

% t0_bar = find(h_T==max(h_T));             % timing phase
t0_bar = length(gm);
r_c_prime = r_c_prime(t0_bar:end);
x = downsample(r_c_prime,4);

r_gm = xcorr(gm,gm);                     % crosscorr ??
rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

% save('Receiver_ab.mat')

%% RECEIVER a
%clc; close all; clear global; clearvars;
% load('rec_input.mat')
% load('Receiver_ab.mat')

M1 = 5;
M2 = 0;
D = 2;
[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

detected = equalization_LE(x, c_opt, M1, D, max(conv(c_opt, h_T)));

[Pe, errors] = SER(in_bits(1:length(detected)), detected);

% %% RECEIVER b
% M1 = 5;
% M2 = 3;
% D = 2;
% [c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
% 