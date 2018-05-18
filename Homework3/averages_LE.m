clc; close all; clear global; clearvars;

load('Useful.mat', 'in_bits', 'qc');

SNR_vect = 8:14;
sigma_a = 2;
M = 4;
gm = conj(qc(end:-1:1));
h = conv(qc,gm);
t0_bar = find(h == max(h));
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
r_gm = xcorr(gm,gm);
realizations = 1:10;
Pe_LE_avg = zeros(length(SNR_vect),1);
Pe_LE = zeros(length(realizations),1);

for i=1:length(SNR_vect)
	Pe_LE = zeros(length(SNR_vect),1);
	for k=1:length(realizations)
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		[r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
		w = wgn(length(r_c),1, 10*log10(sigma_w), 'complex');
		r_c = r_c + w;
		r_c_prime = filter(gm,1,r_c);
		r_c_prime = r_c_prime(t0_bar:end);
		x = downsample(r_c_prime,4);
		rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
		M1 = 7;
		M2 = 0;
		D = 6;
		[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
		detected = equalization_LE(x, c_opt, M1, D, max(conv(c_opt, h_T)));
		[Pe_LE(k),~] = SER(in_bits(1:length(detected)), detected);
	end
	Pe_LE_avg(i) = sum(Pe_LE)/length(Pe_LE);
end

figure();
semilogy(SNR_vect, Pe_LE_avg, 'b--');
grid on;
ylim([10^-4 10^-1]); xlim([8 14]);

% save('PE_LE_avgs.mat', 'Pe_LE_avg');