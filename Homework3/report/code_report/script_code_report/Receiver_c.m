clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex');

% Load input, noise and filter
load('Useful.mat');
load('GAA_filter.mat');

% Channel SNR
snr_db = 10;
snr_lin = 10^(snr_db/10);

sigma_a = 2;	% Input variance

[r_c, sigma_w, ~] = channel_sim(in_bits, snr_db, sigma_a);
r_c = r_c + w(:,3);

r_c_prime = filter(g_AA , 1, r_c);	% Filtering using antialiasing

qg_up = conv(qc, g_AA);
qg_up = qg_up.';

t0_bar = find(qg_up == max(qg_up));				% Timing phase
x = downsample(r_c_prime(t0_bar:end), 2);
qg = downsample(qg_up(1:end), 2);
g_m = conj(flipud(qg));							% Matched filter

x_prime = filter(g_m, 1, x);
x_prime = x_prime(13:end);

h = conv(qg, g_m);

r_g = xcorr(conv(g_AA, g_m));		% AA and MF crosscorrelation
N0 = (sigma_a * E_qc) / (4 * snr_lin);
rw_tilde = N0 .* downsample(r_g, 2);

% Parameters for Equalizer
N1 = floor(length(h)/2);
N2 = N1;
M1 = 5;
D = 4;
M2 = N2 + M1 - 1 - D;

[c, Jmin] = WienerC_frac(h, rw_tilde, sigma_a, M1, M2, D, N1, N2);
psi = conv(h,c);					% Overall impulse response
psi_down = downsample(psi(2:end),2); % The b filter act at T
b = -psi_down(find(psi_down == max(psi_down)) + 1:end); 
x_prime = x_prime/max(psi);			% Normalization
detected = equalization_pointC(x_prime, c_opt, b, D);
detected = detected(1:end-D);
in_bits_2 = in_bits(1:length(detected));
errors = length(find(in_bits_2~=detected(1:length(in_bits_2))));
Pe = errors/length(in_bits_2);