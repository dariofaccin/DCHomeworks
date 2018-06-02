clc; close all; clear global; clearvars;

load('Input_symbols.mat');

symbols = bitmap(x).';
SNR_vect = 0:14;
sigma_a = 2;	% Input variance
Pbit_OFDM_uncode = zeros(length(SNR_vect),1);
M = 512;		% Sub-channels
Npx = 18;		% Prefix length
tic
parfor k=1:length(SNR_vect)
	snr_db = SNR_vect(k);
	snr_lin = 10^(snr_db/10);
	[r_c, sigma_w, g_srrc, g, t0] = channel_OFDM(symbols, snr_db, sigma_a, Npx);
	G = fft(g,512).';
	a_matrix = reshape(r_c(1:end-mod(length(r_c),M+Npx)), M+Npx, []);
	rn = a_matrix(Npx+1:end,:);
	x_k = fft(rn);
	K_i = 1./G;
	y_matrix = x_k.*K_i;
	y = reshape(y_matrix,1,[]);
	hard_d = zeros(length(y),1);
	for i=1:length(hard_d)
		hard_d(i) = QPSK_detector(y(i));
	end
	hard_bits = ibmap(hard_d);
	nerr = length(find(x(1:length(hard_bits))~=hard_bits));
	Pbit_OFDM_uncode(k) = nerr/length(hard_bits);
end
toc