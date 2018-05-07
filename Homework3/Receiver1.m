clc; close all; clear global; clearvars;

in_sig = PNSeq(1023);
in_sig = in_sig(1:length(in_sig)-1);

in_bits = bitmap(in_sig);

sigma_a = 2;
snr_db = 10;
snr_lin = 10^(10/10);

[ch_out, sigma_w] = channel_sim(in_bits, snr_lin, sigma_a);