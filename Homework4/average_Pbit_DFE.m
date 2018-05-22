clc; close all; clear global; clearvars;

load('Useful.mat', 'qc');
load('Input_symbols.mat');

input_bits = x;
SNR_vect = 0:14;
sigma_a = 2;	% Input variance
M = 4;			% Constellation cardinality

gm = conj(qc(end:-1:1));		% Matched filter: complex conjugate of qc
h = conv(qc,gm);				% Impulse response

t0_bar = find(h == max(h));		% Timing phase: peak of h
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);

r_gm = xcorr(gm,gm);			% Matched filter autocorrelation
realizations = 1:1;
Pbit_DFE_avg = zeros(length(SNR_vect),1);
Pe_H_DFE_avg = zeros(length(SNR_vect),1);
Pbit_DFE = zeros(length(realizations),1);
Pe_H_DFE = zeros(length(realizations),1);

N2 = floor(length(h_T)/2);
N1 = N2;
M1 = 5;
D = 4;
M2 = N2 + M1 - 1 - D;
% parpool(2);
tic
max_real = length(realizations);

for i=1:length(SNR_vect)
	Pbit_DFE = zeros(length(realizations),1);
	Pe_H_DFE = zeros(length(realizations),1);
	for k=1:max_real
		snr_db = SNR_vect(i);
		snr_lin = 10^(snr_db/10);
		[r_c, sigma_w, qc] = channel_sim(symbols_ak, snr_db, sigma_a);
		w = wgn(length(r_c),1, 10*log10(sigma_w), 'complex');
		r_c = r_c + w;
		r_c_prime = filter(gm,1,r_c);
		r_c_prime = r_c_prime(t0_bar:end);
		x_aa = downsample(r_c_prime,4);
		rw_tilde = sigma_w/4 .* downsample(r_gm, 4);
		[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
		psi = conv(c_opt, h_T);
		b = - psi(end - M2 + 1:end);
		y_hat = x_aa/max(psi);
		detected = equalization_DFE(y_hat, c_opt, b, M1, M2, D);
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
		Pbit_DFE(k) = nerr/length(decoded);
		Pe_H_DFE(k) = length(find(symbols_ak(1:length(hard_d))~=hard_d))/length(hard_d);

	end
	Pbit_DFE_avg(i) = sum(Pbit_DFE)/length(Pbit_DFE);
	Pe_H_DFE_avg(i) = sum(Pe_H_DFE)/length(Pe_H_DFE);
end
toc

figure();
semilogy(SNR_vect, Pbit_DFE_avg, 'b', 'Marker', '^');
hold on; grid on;
semilogy(SNR_vect, Pe_H_DFE_avg/2, 'r', 'Marker', 'o');
ylim([10^-5 10^-1]); xlim([0 14]);
legend('Coded QPSK + DFE', 'Uncoded QPSK + DFE');

% save('Pbit_DFE_avgs.mat', 'Pbit_DFE_avg', 'Pe_H_DFE_avg');