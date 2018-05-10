clc; close all; clear global; clearvars;

load('Useful.mat', 'in_bits', 'qc');

SNR_vect = [8 11 14];

sigma_a = 2;
M = 4;

gm = conj(qc(end:-1:1));
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);
h_T = downsample(h,4);
t0_bar = length(gm);
r_gm = xcorr(gm,gm);

N2 = floor(length(h_T)/2);

printmsg_delete = '';
JminDFE = zeros(length(SNR_vect), 30, 30);

for snr=1:length(SNR_vect)
    channel_snr_db = SNR_vect(snr);
    for M1=1:30
        for D=1:30
            printmsg = sprintf('snr = %d, M1 = %d, D = %d\n', channel_snr_db, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            [~, sigma_w, ~] = channel_sim(in_bits, channel_snr_db, sigma_a);
            M2 = N2 + M1 - 1 - D;
            r_w = sigma_w/4 .* downsample(r_gm, 4);
            [~, JminDFE(snr,M1,D)] = Adaptive_DFE(h_T, r_w, sigma_a, M1, M2, D);
        end
    end
end

for i = 1:length(SNR_vect)
    figure, mesh(1:30, 1:30, reshape(abs(JminDFE(i, :, :)), size(JminDFE(i, :, :), 2), size(JminDFE(i, :, :), 3)))
    title(strcat('Jmin for LE, snr= ', num2str(SNR_vect(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin [dB]')
end

[min, idx] = min(JminDFE(:));

[idx_snr, idx_m1, idx_d] = ind2sub(size(JminDFE), idx);

save('JminDFE.mat', 'JminDFE', 'idx_m1', 'idx_d', 'M2');