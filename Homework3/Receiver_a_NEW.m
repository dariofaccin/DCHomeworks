clc; close all; clear global; clearvars;

% Load input and noise
load('Useful.mat');

% Channel SNR
snr_db = 10;
snr_lin = 10^(snr_db/10);

sigma_a = 2;

% Channel: NOISE IS ADDED AFTERWARDS
[r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);

r_c = r_c + w(:,3);

% Matched filter
gm = conj(qc(end:-1:1));

% Impulse response
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);

% Downsampling impulse response
h_T = downsample(h,4);

% Filtering received signal
r_c_prime = filter(gm,1,r_c);

% Determining timing phase
t0_bar = length(gm);

% Remove "transient" and downsample received signal
r_c_prime = r_c_prime(t0_bar:end);
x = downsample(r_c_prime,4);

% Filter autocorrelation
r_gm = xcorr(gm,gm);
rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

% Parameters for Linear Equalizer
M1 = 5;
M2 = 0;
D = 2;
[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

psi = conv(c_opt, h_T);

figure
subplot(121), stem(0:length(c_opt)-1,abs(c_opt)), hold on, grid on
title('|c|'), xlabel('n');
subplot(122), stem(0:length(psi)-1,abs(psi)), grid on
title('|\varpsi|'), xlabel('n');

detected = equalization_LE(x, c_opt, M1, D, max(psi));

nerr = length(find(in_bits(1:length(detected))~=detected));
Pe = nerr/length(in_bits(1:length(detected)));


% [Pe, errors] = SER(in_bits(1:length(detected)), detected);