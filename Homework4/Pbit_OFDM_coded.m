clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

input_bits = x;
SNR_vect = 10:14;
numb = length(SNR_vect);
fprintf('Number of iterations = %d\n', numb );
sigma_a = 2;	% Input variance

realizations = 1:1;
Pbit_OFDM_avg = zeros(length(SNR_vect),1);
% Pe_H_OFDM_avg = zeros(length(SNR_vect),1);
Pbit_OFDM = zeros(length(realizations),1);
% Pe_H_OFDM = zeros(length(realizations),1);
t0 = 21;
M = 512;
Npx = 8;

for i=1:length(SNR_vect)
	Pbit_OFDM = zeros(length(realizations),1);
	for k=1:length(realizations)
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		[r_c, sigma_w, g_srrc, g, t0] = channel_OFDM(symbols_ak, snr_db, sigma_a);
        G = fft(g,512).';
        a_matrix = reshape(r_c(1:end-mod(length(r_c),M+Npx)), M+Npx, []);
        rn = a_matrix(Npx+1:end,:);
        x_k = fft(rn);
        K_i = 1./G;
        y_matrix = x_k.*K_i;
        sigma_i = 0.5*sigma_w*M*abs(K_i).^2;
        llr_real = -2*real(y_matrix).*sigma_i.^(-1);
        llr_imag = -2*imag(y_matrix).*sigma_i.^(-1);
        llr_real_ar = reshape(llr_real, [], 1);
        llr_imag_ar = reshape(llr_imag, [], 1);
        llr = zeros(numel(llr_real) + numel(llr_imag), 1);
        llr(1:2:end) = llr_real_ar;
        llr(2:2:end) = llr_imag_ar;
        llr = deinterleaver(llr);
        tic
        dec_bits = LDPC_decoder(llr).';
        toc
        nerr = length(find(x(1:length(dec_bits))~=dec_bits));
        Pbit_OFDM(k) = nerr/length(dec_bits);
	end
	Pbit_OFDM_avg(i) = sum(Pbit_OFDM)/length(Pbit_OFDM);
end
toc

figure();
semilogy(SNR_vect, Pbit_OFDM_avg, 'b', 'Marker', '^');
hold on; grid on;
ylim([10^-5 10^-1]); xlim([SNR_vect(1) SNR_vect(end)]);
legend('Coded OFDM');

save('Pbit_OFDM_avgs_coded','Pbit_OFDM_avg')