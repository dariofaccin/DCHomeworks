clc; close all; clear global; clearvars;

% Initial bits
L = 2^15-1;
x = PNSeq(L);
N = floor(length(x)/32400);        % to be encoded by H = 32400x64800
x = x(1:N*32400);

% Matlab LDPC Encoder
h = comm.LDPCEncoder;
encoded = LDPC_encoder(x,h,N);

interleaved = interleaver(encoded);
symbols_ak = bitmap(interleaved.').';

% save('Input_symbols.mat', 'symbols_ak', 'x');