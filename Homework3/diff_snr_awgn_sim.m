clc; close all; clear global; clearvars;

load('Useful.mat');

sigma_a = 2;
M = 4;
SNR_vect = 8:14;
Pe_AWGN_sim = zeros(length(SNR_vect),1);
errors = zeros(length(in_bits),1);

for i=1:length(SNR_vect)
	snr_db = SNR_vect(i);
	sr_lin = 10^(snr_db/10);
	r_c = in_bits + w(1:length(in_bits),i);
	detected = zeros(length(r_c),1);
	for k=1:length(in_bits)
		detected(k) = QPSK_detector(r_c(k));
	end
	errors(i) = length(find(in_bits(1:length(detected))~=detected));
    Pe_AWGN_sim(i) = errors(i)/length(in_bits(1:length(detected)));
end

figure, semilogy(SNR_vect, Pe_AWGN_sim);

save('Pe_AWGN_sim.mat', 'Pe_AWGN_sim');