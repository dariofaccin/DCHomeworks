clc; close all; clear global; clearvars;

% Input
load('Useful.mat');
SNR_vect = 8:14;
Pe_FBA = zeros(length(SNR_vect),1);
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

printmsg_delete = '';
parpool('local',2);

for i=1:length(SNR_vect)
    printmsg = sprintf('snr = %d\n', SNR_vect(i));
    fprintf([printmsg_delete, printmsg]);
    printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
    
    snr_db = SNR_vect(i);
    snr_lin = 10^(snr_db/10);
    [r_c, sigma_w, qc] = channel_sim(in_bits, snr_db, sigma_a);
    r_c = r_c + w(:,i);
    
    % Filtering received signal
    r_c_prime = filter(gm,1,r_c);
    
    r_c_prime = r_c_prime(t0_bar:end);
    x = downsample(r_c_prime,4);

    rw_tilde = sigma_w/4 .* downsample(r_gm, 4);

    M1 = 5;
    M2 = 0;
    D = 1;
    
    [c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);
    
    psi = conv(c_opt, h_T);
    psi = psi/max(psi);
    
    y = conv(x, c_opt);
    y = y/max(psi);
    indexD = find(psi == max(psi));
    L1 = 2; L2 = 2;
    
    detected = FBA(y, psi(indexD-L1:indexD+L2), L1, L2);

    errors(i) = length(find(in_bits(1:length(detected))~=detected));
    Pe_FBA(i) = errors(i)/length(in_bits(1:length(detected)));
end

figure, semilogy(SNR_vect, Pe_FBA, SNR_vect, awgn_bound,'r--');

save('Pe_FBA.mat','Pe_FBA');