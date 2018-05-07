clc; close all; clear global; clearvars;

in_sig = PNSeq(1023);
in_sig = in_sig(1:length(in_sig)-1);

in_bits = bitmap(in_sig);

sigma_a = 2;
snr_db = 10;
snr_lin = 10^(10/10);

[ch_out, sigma_w, qc] = channel_sim(in_bits, snr_lin, sigma_a);

q_mf = qc(end:-1:1);

y = filter(q_mf,1,ch_out);

y = downsample(y,2);

detected = zeros(length(y),1);

for i=1:length(detected)
    detected(i) = QPSK_detector(y(i));
end

numerr = 0;
for i=1:length(detected)
    if ( detected(i) ~= in_bits(i))
        numerr = numerr +1;
    end
end
