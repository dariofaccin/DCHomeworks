clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
SNR_vect = [8 9 10 11 12 13 14];
Pe_AA_NOGM = zeros(length(SNR_vect),1);
errors = zeros(length(SNR_vect),1);
sigma_a = 2;

%% AA filter

Fpass = 0.2;             % Passband Frequency
Fstop = 0.3;             % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.01;            % Stopband Attenuation
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);

% Calculate the coefficients using the FIRPM function.
g_AA  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(g_AA);

qg_up = conv(qc, g_AA);
qg_up = qg_up.';
t0_bar = find(qg_up == max(qg_up));

qg = downsample(qg_up, 2);

h = qg;

r_g = xcorr(g_AA);
N1 = floor(length(h)/2);
N2 = 12;
M1 = 9;
D = 4;
M2 = N2 + M1 - 1 - D;

for i=1:length(SNR_vect)
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, ~] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);

	r_c_prime = filter(g_AA , 1, r_c);

	x = downsample(r_c_prime(t0_bar:end), 2);

	x_prime = x;
	
	%% Equalization and symbol detection

	N0 = (sigma_a * 1) / (4 * snr_lin);
	r_w = N0 * downsample(r_g, 2);

	[c, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2);
	psi = conv(h,c);

	psi_down = downsample(psi(2:end),2); % The b filter act at T
	b = -psi_down(find(psi_down == max(psi_down)) + 1:end); 

	detected = equalization_pointC(x_prime, c, b, D);
	
    [Pe_AA_NOGM(i), errors(i)] = SER(in_bits(1:length(detected)), detected(1:end));
end

figure();
semilogy(SNR_vect, Pe_AA_NOGM, 'k'); grid on;
ylim([10^-4 10^-1]); xlim([8 14]);
legend('AAF+DFE@$\frac{T}{2}$'); set(legend,'Interpreter','latex');

% save('Pe_AA_NOGM.mat','Pe_AA_NOGM');