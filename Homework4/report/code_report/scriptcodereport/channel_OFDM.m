function[output, sigma_w, g_srrc, tot_ds, t0, tot] = channel_OFDM(input, snr, sigma_a, Npx)

snr_db = snr;
snr_lin = 10^(snr_db/10);
M = 512;

a_pad = [input; ones(M - mod(length(input), M), 1) * (1+1i)];
a_matrix = reshape(a_pad, M, []); % Should mantain columnwise order

A_matrix = ifft(a_matrix);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix];

r = reshape(A_matrix, [], 1);

Q = 4;
in_upsampled = upsample(r, Q);

% Square-root raised cosine
ro = 0.0625;
span = 30;
sps = 4;
g_srrc = rcosdesign(ro, span, sps, 'sqrt');

in_after_srrc = filter(g_srrc,1,in_upsampled);

alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
qc = [0; 0; 0; 0; 0; qc(qc >=max(qc)*10^(-2))];

in_after_qc = filter(qc,1,in_after_srrc);
all_ch = conv(g_srrc, conv(g_srrc,qc));
E_tot = sum(conv(g_srrc,qc).^2);
sigma_w = sigma_a/M * E_tot / snr_lin;
in_after_qc = in_after_qc + wgn(length(in_after_qc),1,10*log10(sigma_w),'complex');

tot = all_ch(abs(all_ch)>=(max(abs(all_ch))*1e-2));
tot_ds = downsample(tot, 4);

in_after_srrc = filter(g_srrc, 1, in_after_qc);
t0 = find(in_after_srrc==max(in_after_srrc));

in_after_srrc = in_after_srrc(t0:end);
in_after_srrc = downsample(in_after_srrc,4);

output = in_after_srrc;
end