function[output, sigma_w, g_srrc, tot_ds, t0] = channel_OFDM(input, snr, sigma_a)

snr_db = snr;
snr_lin = 10^(snr_db/10);

M = 512;
Npx = 5;

a_pad = [input; ones(M - mod(length(input), M), 1) * (1+1i)];
a_matrix = reshape(a_pad, M, []); % it should mantain columnwise order

A_matrix = ifft(a_matrix);
A_matrix = [A_matrix(M-Npx+1:M, :); A_matrix]; % very powerful operationt

r = reshape(A_matrix, [], 1);

Q = 4;
in_upsampled = upsample(r, Q);

% Square-root raised cosine
% N    = 26;         % Order
% Fc   = 0.46;       % Cutoff Frequency
% TM   = 'Rolloff';  % Transition Mode
% R    = 0.0625;     % Rolloff
% DT   = 'sqrt';     % Design Type
% Beta = 0.7;        % Window Parameter
% win = kaiser(N+1, Beta);

% g_srrc  = firrcos(N, Fc, R, 2, TM, DT, [], win);

g_srrc = rcosdesign(0.0625, 12, 4, 'sqrt');

in_after_srrc = filter(g_srrc,1,in_upsampled);

alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
qc = [0; 0; 0; 0; 0; qc(qc >=max(qc)*10^(-2))];

in_after_qc = filter(qc,1,in_after_srrc);
all_ch = conv(g_srrc, conv(g_srrc,qc));
E_tot = sum(all_ch.^2);
sigma_w = sigma_a/M * E_tot / snr_lin;
tot = all_ch(3+25-1:end-25);
tot_ds = downsample(tot, 4);

t0 = 29;
in_after_srrc = filter(g_srrc, 1, in_after_qc);

in_after_srrc = in_after_srrc(t0:end);
in_after_srrc = downsample(in_after_srrc,4);

output = in_after_srrc;