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

printmsg_delete = '';
JminLE = zeros(length(SNR_vect), 30, 30);

for snr=1:length(SNR_vect)
    channel_snr_db = SNR_vect(snr);
    for M1=1:30
        for D=1:30
            printmsg = sprintf('snr = %d, M1 = %d, D = %d\n', channel_snr_db, M1, D);
            fprintf([printmsg_delete, printmsg]);
            printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
            [~, sigma_w, ~] = channel_sim(in_bits, channel_snr_db, sigma_a);
            M2 = 0;
            r_w = sigma_w/4 .* downsample(r_gm, 4);
            [~, JminLE(snr,M1,D)] = Adaptive_DFE(h_T, r_w, sigma_a, M1, M2, D);
        end
    end
end

for i = 1:length(SNR_vect)
    figure, mesh(1:30, 1:30, reshape(abs(JminLE(i, :, :)), size(JminLE(i, :, :), 2), size(JminLE(i, :, :), 3)))
    title(strcat('Jmin for LE, snr= ', num2str(SNR_vect(i))))
    xlabel('D'), ylabel('M1'), zlabel('Jmin [dB]')
end

[min, idx] = min(JminLE(:));

[idx_snr, idx_m1, idx_d] = ind2sub(size(JminLE), idx);

save('JminLE.mat', 'JminLE', 'idx_m1', 'idx_d');