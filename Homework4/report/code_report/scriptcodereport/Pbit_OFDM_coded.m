clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

SNR_vect = 0.8:0.05:1.4;
sigma_a = 2;	% Input variance

Pbit_OFDM_code = zeros(length(SNR_vect),1);
M = 512;		% Sub-channels
Npx = 18;		% Prefix length
tic
parfor i=1:length(SNR_vect)
	snr_db = SNR_vect(i);
	snr_lin = 10^(snr_db/10);
    % simulate transmission over the OFDM channel
	[r_c, sigma_w, g_srrc, g, t0] = channel_OFDM(symbols_ak, snr_db, sigma_a, Npx);
	G = fft(g,512).';
	a_matrix = reshape(r_c(1:end-mod(length(r_c),M+Npx)), M+Npx, []);
	rn = a_matrix(Npx+1:end,:);
    % DFT
	x_k = fft(rn);
	K_i = 1./G;
    % output of the channel
	y_matrix = x_k.*K_i;
	sigma_i = 0.5*sigma_w*M*abs(K_i).^2;
    % LLR
	llr_real = -2*real(y_matrix).*sigma_i.^(-1);
	llr_imag = -2*imag(y_matrix).*sigma_i.^(-1);
	llr_real_ar = reshape(llr_real, [], 1);
	llr_imag_ar = reshape(llr_imag, [], 1);
	llr = zeros(numel(llr_real) + numel(llr_imag), 1);
	llr(1:2:end) = llr_real_ar;
	llr(2:2:end) = llr_imag_ar;
	llr = deinterleaver(llr);
    % LDPC decoder
	dec_bits = LDPC_decoder(llr).';
	nerr = length(find(x(1:length(dec_bits))~=dec_bits));
	Pbit_OFDM_code(i) = nerr/length(dec_bits);
end
toc