clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')

% Load input and noise
load('Useful.mat');

% Channel SNR
snr_db = 10;
snr_lin = 10^(snr_db/10);

sigma_a = 2;	% Input variance

[r_c, sigma_w, ~] = channel_sim(in_bits, snr_db, sigma_a);
r_c = r_c + w(:,3);

gm = conj(qc(end:-1:1));	% Matched filter: complex conjugate of qc

figure()
stem(abs(gm));
xlabel('$m\frac{T}{4}$');
ylabel('$g_m$')
xlim([1 length(gm)]);
grid on

h = conv(qc,gm);			% Impulse response
h = h(h>max(h)/100);
h = h(3:end-2);

h_T = downsample(h,4);		% Downsampling impulse response
r_c_prime = filter(gm,1,r_c);	% Filtering received signal
t0_bar = find(h == max(h));		% Determining timing phase

r_c_prime = r_c_prime(t0_bar:end);	% Remove "transient"
x = downsample(r_c_prime,4);		% Downsample received signal

r_gm = xcorr(gm,gm);			% Filter autocorrelation
rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

% Parameters for Linear Equalizer
M1 = 7;
M2 = 0;
D = 6;
[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

psi = conv(c_opt, h_T);	% Overall impulse response

figure
subplot(121), stem(0:length(c_opt)-1,abs(c_opt)), hold on, grid on
ylabel('$|c|$'), xlabel('n'); xlim([0 7]);
subplot(122), stem(-2:length(psi)-3,abs(psi)), grid on
ylabel('$|\psi|$'), xlabel('n'); xlim([-2 8]);

detected = equalization_LE(x, c_opt, M1, D, max(psi));

nerr = length(find(in_bits(1:length(detected))~=detected));
Pe = nerr/length(detected);