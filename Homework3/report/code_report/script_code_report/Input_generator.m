clc; close all; clear global; clearvars;

% CREATE INPUT BITS, QC AND NOISE

alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
qc = [0; 0; 0; 0; 0; qc(qc >=max(qc)/100)];
E_qc = sum(qc.^2);

length_seq = 2^20-1;

SNR_vect = 8:14;
SNR_lin = 10.^(SNR_vect ./ 10);

sigma_w = zeros(length(SNR_vect),1);
sigma_a = 2;

w = zeros(2*length_seq-2, 7);

for i=1:length(SNR_vect)
    sigma_w(i) = E_qc * sigma_a / SNR_lin(i);
    w(:,i) = wgn(2*length_seq-2,1, 10*log10(sigma_w(i)), 'complex');
end

in_seq = PNSeq(length_seq);
in_seq = in_seq(1:end-1);
in_bits = bitmap(in_seq);
% save('Useful.mat', 'w', 'in_bits', 'qc');