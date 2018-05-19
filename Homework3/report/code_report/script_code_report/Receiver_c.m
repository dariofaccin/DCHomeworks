clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')    % latex format

% Load input, noise and filter
load('Useful.mat');
load('GAA_filter.mat');

% Channel SNR
snr_db = 10;
snr_lin = 10^(snr_db/10);

sigma_a = 2;

% Channel: NOISE IS ADDED AFTERWARDS
[r_c, sigma_w, ~] = channel_sim(in_bits, snr_db, sigma_a);
s_c = r_c;                  % Useful noise
r_c = r_c + w(:,3);

r_c_prime = filter(g_AA , 1, r_c);

qg_up = conv(qc, g_AA);
qg_up = qg_up.';
%freqz(qg_up, 1,'whole');
figure, stem(qg_up), title('convolution of $g_AA$ and $q_c$'), xlabel('nT/4')

%% Timing phase and decimation

t0_bar = find(qg_up == max(qg_up));
x = downsample(r_c_prime(t0_bar:end), 2);

qg = downsample(qg_up(1:end), 2);
g_m = conj(flipud(qg));

[GM, ff] =  freqz(g_m,1,'whole');
figure, plot(2*ff/(pi),20*log10(abs(GM))), xlim([0 2]),
ylabel('$|G_M|$ [dB]')
ylim([-40 10]);
xlabel('f/T')
grid on;

figure, stem(g_m), title('$g_m$'), xlabel('nT/2')

x_prime = filter(g_m, 1, x);
x_prime = x_prime(13:end);

h = conv(qg, g_m);
h = h(h ~= 0);

%% Equalization and symbol detection

r_g = xcorr(conv(g_AA, g_m));
N0 = (sigma_a * 1) / (4 * snr_lin);
r_w = N0 * downsample(r_g, 2);

figure, stem(r_w), title('$r_w$'), xlabel('nT/2')
figure, stem(r_g), title('$r_g$'), xlabel('nT/2')

N1 = floor(length(h)/2);
N2 = N1;

M1 = 5;
D = 4;
M2 = N2 + M1 - 1 - D;

[c, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2);
psi = conv(h,c);

figure, stem(c), title('c'), xlabel('nT/2')
figure, stem(abs(psi)), title('|$\psi$|'), xlabel('nT/2')

psi_down = downsample(psi(2:end),2); % The b filter act at T
b = -psi_down(find(psi_down == max(psi_down)) + 1:end); 

figure, stem(b), title('b'), xlabel('nT')
detected = equalization_pointC(x_prime, c, b, D);

nerr = length(find(in_bits(D:length(detected)+D-1)~=detected));
Pe = nerr/length(detected);