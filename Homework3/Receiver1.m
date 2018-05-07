clc; close all; clear global; clearvars;

% Input signal
load('rec_input.mat');

snr_db = 10;
snr_lin = 10^(snr_db/10);

% Matched filter
q_mf = qc(end:-1:1);

% Filtering through C
y = filter(q_mf,1,r_c);
y = downsample(y,4);

detected = zeros(length(y),1);

for i=1:length(detected)
    detected(i) = QPSK_detector(y(i));
end

numerrs = 0;
for i=1:length(detected)
    if ( detected(i) ~= in_bits(i) )
        numerrs = numerrs + 1;
    end
end