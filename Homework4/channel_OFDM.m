function [output, sigma_w, g_srrc, all_ch] = channel_OFDM(input, snr, sigma_a)

M = 512;

snr_db = snr;
snr_lin = 10^(snr_db/10);

Q = 4;
in_upsampled = upsample(input, Q);

% Square-root raised cosine
N    = 26;         % Order
Fc   = 0.46;       % Cutoff Frequency
TM   = 'Rolloff';  % Transition Mode
R    = 0.0625;     % Rolloff
DT   = 'sqrt';     % Design Type
Beta = 0.5;        % Window Parameter
win = kaiser(N+1, Beta);

g_srrc  = firrcos(N, Fc, R, 2, TM, DT, [], win);

in_after_srrc = filter(g_srrc,1,in_upsampled);

alpha = 0.67;
beta = 0.7424;
qc_num = [0 0 0 0 0 beta];
qc_denom = [1 -alpha];
qc = impz(qc_num, qc_denom);
qc = [0; 0; 0; 0; 0; qc(qc >=max(qc)*10^(-2))];

in_after_qc = filter(qc,1,in_after_srrc);

all_ch = conv(g_srrc, qc);

E_tot = sum(all_ch.^2);

sigma_w = sigma_a/M * E_tot / snr_lin;

output = in_after_qc;