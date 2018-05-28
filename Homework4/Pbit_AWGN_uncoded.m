clc; close all; clear global; clearvars;

load('Input_symbols.mat');

symbols = bitmap(x).';

SNR_vect = 0:14;
sigma_a = 2;

Pbit_AWGN_uncode = zeros(length(SNR_vect),1);

tic
parfor i=1:length(SNR_vect)
	snr_db = SNR_vect(i);
	snr_lin = 10^(snr_db/10);
	sigma_w = sigma_a / snr_lin;
	w = wgn(length(symbols),1, 10*log10(sigma_w), 'complex');
	r_c = symbols + w;
	r_c_prime = r_c;
	detected = r_c_prime;
	hard_d = zeros(length(detected),1);
	for u=1:length(detected)
		hard_d(u) = QPSK_detector(detected(u));
	end
	hard_bits = ibmap(hard_d);
	Pbit_AWGN_uncode(i) = length(find(x(1:length(hard_bits))~=hard_bits))/length(hard_bits);
end
toc

figure();
semilogy(SNR_vect, Pbit_AWGN_uncode, 'b', 'Marker', '^');
grid on;
ylim([10^-5 10^-1]); xlim([0 14]);
legend('Uncoded AWGN');

% save('Pbit_AWGN_uncoded.mat', 'Pbit_AWGN_uncode');