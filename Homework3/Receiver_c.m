clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')

% Load input, noise and filter
load('Useful.mat');
load('GAA_filter.mat');

% Anti aliasing filter
[G_AA, f] =  freqz(g_AA,1,'whole');
figure, plot(4*f/(2*pi),20*log10(abs(G_AA))), xlim([0 2]),
ylabel('$|G_{AA}|$ [dB]')
ylim([-40 5]);
xlabel('f/T')
grid on;

% Channel SNR
snr_db = 10;
snr_lin = 10^(snr_db/10);

sigma_a = 2;	% Input variance

[r_c, sigma_w, ~] = channel_sim(in_bits, snr_db, sigma_a);
r_c = r_c + w(:,3);

r_c_prime = filter(g_AA , 1, r_c);	% Filtering using antialiasing

qg_up = conv(qc, g_AA);
qg_up = qg_up.';

% figure, stem(qg_up), title('convolution of $g_AA$ and $q_c$'), xlabel('nT/4')

%% Timing phase and decimation

t0_bar = find(qg_up == max(qg_up));
x = downsample(r_c_prime(t0_bar:end), 2);

qg = downsample(qg_up(1:end), 2);
g_m = conj(flipud(qg));

[GM, ff] =  freqz(g_m,1,'whole');
figure, plot(2*ff/(2*pi),20*log10(abs(GM)))
xlim([0 2])
ylabel('$|G_M|$ [dB]')
ylim([-15 10]);
xlabel('f/T')
grid on;

% figure, stem(g_m), title('$g_m$'), xlabel('nT/2')

x_prime = filter(g_m, 1, x);
x_prime = x_prime(13:end);

h = conv(qg, g_m);
h = h(h ~= 0);

%% Equalization and symbol detection

r_g = xcorr(conv(g_AA, g_m));
N0 = (sigma_a * 1) / (4 * snr_lin);
r_w = N0 * downsample(r_g, 2);

% figure, stem(r_w), title('$r_w$'), xlabel('nT/2')
% figure, stem(r_g), title('$r_g$'), xlabel('nT/2')

N1 = floor(length(h)/2);
N2 = N1;
M1 = 10;
D = 4;
M2 = N2 + M1 - 1 - D;

[c_opt, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2);
psi = conv(h,c_opt);

psi_down = downsample(psi(2:end),2); % The b filter act at T
b = -psi_down(find(psi_down == max(psi_down)) + 1:end); 
x_prime = x_prime/max(psi);
detected = equalization_pointC(x_prime, c_opt, b, D);
detected = detected(1:end-D);
in_bits_2 = in_bits(1:length(detected));
errors = length(find(in_bits_2~=detected(1:length(in_bits_2))));
Pe = errors/length(in_bits_2);

%% C
figure, stem(0:length(c_opt)-1,abs(c_opt)), hold on, grid on
ylabel('$|c|$'), xlabel('$n\frac{T}{2}$'); xlim([0 length(c_opt)-1]);
%% B
figure, stem(0:length(b)-1,abs(b)), hold on, grid on
ylabel('$|b|$'), xlabel('n'); xlim([0 length(b)-1]);
%% PSI
figure
stem(-(find(psi==max(psi))-D)+1:length(psi)-(find(psi==max(psi))-D+1)+1,abs(psi))
xlim([-(find(psi==max(psi))-D)+1 length(psi)-(find(psi==max(psi))-D+1)+1]), grid on
ylabel('$|\psi|$'), xlabel('$n\frac{T}{2}$'); 
xlim([-10 15]);