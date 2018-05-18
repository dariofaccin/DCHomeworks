clc; close all; clear global; clearvars;

load('Useful.mat', 'in_bits', 'qc');

SNR_vect = 8:14;
sigma_a = 2;
M = 4;
realizations = 1:10;
Pe_AWGN_SIM_avg = zeros(length(SNR_vect),1);
Pe_AWGN_SIM = zeros(length(realizations),1);
awgn_bound = zeros(length(SNR_vect),1);

for i=1:length(SNR_vect)
	Pe_AWGN_SIM = zeros(length(SNR_vect),1);
	for k=1:length(realizations)
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		r_c = in_bits;
		w = wgn(length(r_c),1, 10*log10(sigma_a / snr_lin), 'complex');
		r_c = r_c + w;
		detected = zeros(length(r_c),1);
		for l=1:length(in_bits)
		detected(l) = QPSK_detector(r_c(l));
		end
		[Pe_AWGN_SIM(k),~] = SER(in_bits(1:length(detected)), detected);
	end
	Pe_AWGN_SIM_avg(i) = sum(Pe_AWGN_SIM)/length(Pe_AWGN_SIM);
	awgn_bound(i) = 4*(1-1/sqrt(M))*qfunc(sqrt(snr_lin/(sigma_a/2)));
end

figure();
semilogy(SNR_vect, Pe_AWGN_SIM_avg, 'g--');
hold on; grid on;
semilogy(SNR_vect, awgn_bound, 'g');
ylim([10^-4 10^-1]); xlim([8 14]);

% save('Pe_AWGN_SIM_avgs.mat', 'Pe_AWGN_SIM_avg', 'awgn_bound');