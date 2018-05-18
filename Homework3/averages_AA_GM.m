clc; close all; clear global; clearvars;

load('Useful.mat', 'in_bits', 'qc');
SNR_vect = 8:14;
sigma_a = 2;
Fpass = 0.2;             % Passband Frequency
Fstop = 0.3;             % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.01;            % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);
g_AA  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(g_AA);
qg_up = conv(qc, g_AA);
qg_up = qg_up.';
t0_bar = find(qg_up == max(qg_up));
qg = downsample(qg_up(1:end), 2);
g_m = conj(flipud(qg));
h = conv(qg, g_m);
h = h(h ~= 0);
r_g = xcorr(conv(g_AA, g_m));
realizations = 1:10;
Pe_AA_GM_avg = zeros(length(SNR_vect),1);
Pe_AA_GM = zeros(length(realizations),1);

for i=1:length(SNR_vect)
	Pe_AA_GM = zeros(length(SNR_vect),1);
	for k=1:length(realizations)
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		[r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
		w = wgn(length(r_c),1, 10*log10(sigma_w), 'complex');
		r_c = r_c + w;
		r_c_prime = filter(g_AA , 1, r_c);
		x = downsample(r_c_prime(t0_bar:end), 2);
		x_prime = filter(g_m, 1, x);
		x_prime = x_prime(13:end);
		r_w = sigma_w/4 .* downsample(r_g, 2);
		N1 = floor(length(h)/2);
		N2 = N1;
		M1 = 5;
		D = 4;
		M2 = N2 + M1 - 1 - D;
		[c, Jmin] = WienerC_frac(h, r_w, sigma_a, M1, M2, D, N1, N2);
		psi = conv(h,c);
		psi_down = downsample(psi(2:end),2); % The b filter act at T
		b = -psi_down(find(psi_down == max(psi_down)) + 1:end); 
		detected = equalization_pointC(x_prime, c, b, D);
		[Pe_AA_GM(k),~] = SER(in_bits(4:length(detected)), detected);
	end
	Pe_AA_GM_avg(i) = sum(Pe_AA_GM)/length(Pe_AA_GM);
end

figure();
semilogy(SNR_vect, Pe_AA_GM_avg, 'k--');
grid on;
ylim([10^-4 10^-1]); xlim([8 14]);

% save('PE_AA_GM_avgs.mat', 'Pe_AA_GM_avg');