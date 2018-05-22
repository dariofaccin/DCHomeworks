clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

input_bits = x;
SNR_vect = 0:14;
sigma_a = 2;	% Input variance
M = 4;			% Constellation cardinality

realizations = 1:1;
Pbit_AWGN_avg = zeros(length(SNR_vect),1);
Pe_H_AWGN_avg = zeros(length(SNR_vect),1);
Pbit_AWGN = zeros(length(realizations),1);
Pe_H_AWGN = zeros(length(realizations),1);

tic
max_real = length(realizations);

for i=1:length(SNR_vect)
	Pbit_AWGN = zeros(length(realizations),1);
	Pe_H_AWGN = zeros(length(realizations),1);
	for k=1:max_real
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		sigma_w = sigma_a / snr_lin;
		w = wgn(length(symbols_ak),1, 10*log10(sigma_w), 'complex');
		r_c = symbols_ak + w;
		r_c_prime = r_c;
		detected = r_c_prime;
        llr = zeros(2*length(detected),1);
        llr(1:2:end) = -2*real(detected)/(sigma_w/2);
        llr(2:2:end) = -2*imag(detected)/(sigma_w/2);
        % Deinterleave and decode bits
        llr_deinter = deinterleaver(llr);
        decoded = LDPC_decoder(llr_deinter);
		hard_d = zeros(length(detected),1);
		for u=1:length(detected)
			hard_d(u) = QPSK_detector(detected(u));
		end
		nerr = length(find(input_bits(1:length(decoded))~=decoded.'));
		Pbit_AWGN(k) = nerr/length(decoded);
		Pe_H_AWGN(k) = length(find(symbols_ak(1:length(hard_d))~=hard_d))/length(hard_d);

	end
	Pbit_AWGN_avg(i) = sum(Pbit_AWGN)/length(Pbit_AWGN);
	Pe_H_AWGN_avg(i) = sum(Pe_H_AWGN)/length(Pe_H_AWGN);
end
toc

figure();
semilogy(SNR_vect, Pbit_AWGN_avg, 'b', 'Marker', '^');
hold on; grid on;
semilogy(SNR_vect, Pe_H_AWGN_avg, 'r', 'Marker', 'o');
ylim([10^-5 10^-1]); xlim([0 14]);
legend('Coded AWGN', 'Uncoded AWGN');

% save('Pbit_AWGN_avgs.mat', 'Pbit_AWGN_avg', 'Pe_H_AWGN_avg');