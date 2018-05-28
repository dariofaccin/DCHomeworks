clc; close all; clear global; clearvars;

load('Input_symbols.mat');

symbols = symbols_ak;

SNR_vect = 2:0.1:3;
sigma_a = 2;

Pbit_AWGN_code = zeros(length(SNR_vect),1);

tic
for i=1:length(SNR_vect)
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
	llr_deinter = deinterleaver(llr);
	decoded = LDPC_decoder(llr_deinter).';
	Pbit_AWGN_code(i) = length(find(x(1:length(decoded))~=decoded))/length(decoded);
end
toc

figure();
semilogy(SNR_vect, Pbit_AWGN_code, 'b', 'Marker', '^');
grid on;
ylim([10^-5 10^-1]); xlim([2 3]);
legend('Coded AWGN');

save('Pbit_AWGN_coded.mat', 'Pbit_AWGN_code');