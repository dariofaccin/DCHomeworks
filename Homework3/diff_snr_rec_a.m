clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
% load('JminLE.mat');
SNR_vect = 8:14;
Pe_LE = zeros(length(SNR_vect),1);
errors = zeros(length(SNR_vect),1);
awgn_bound = zeros(length(SNR_vect),1);

sigma_a = 2;
M = 4;

gm = conj(qc(end:-1:1));
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
t0_bar = length(gm);
r_gm = xcorr(gm,gm);

for i=1:length(SNR_vect)
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);
    
    % Filtering received signal
    r_c_prime = filter(gm,1,r_c);
    
    r_c_prime = r_c_prime(t0_bar:end);
    x = downsample(r_c_prime,4);

    rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

    M1 = 7;
    M2 = 0;
    D = 6;
    [c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

    detected = equalization_LE(x, c_opt, M1, D, max(conv(c_opt, h_T)));

    errors(i) = length(find(in_bits(1:length(detected))~=detected));
    Pe_LE(i) = errors(i)/length(in_bits(1:length(detected)));
    
    awgn_bound(i) = 4*(1-1/sqrt(M))*qfunc(sqrt(snr_lin/(sigma_a/2)));
end


figure();
semilogy(SNR_vect, Pe_LE, 'b--'); hold on; grid on;
ylim([10^-4 10^-1]); xlim([8 14]);
semilogy(SNR_vect, awgn_bound, 'g');

save('Pe_LE.mat','Pe_LE', 'awgn_bound');