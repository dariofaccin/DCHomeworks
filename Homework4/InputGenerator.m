clc; close all; clear global; clearvars;

% Initial bits
L = 2^20-1;
x = PNSeq(L);
% Matlab LDPC Encoder
H = comm.LDPCEncoder();

sstep = 32400;
numbits = floor(length(x) / sstep) * sstep;
x = x(1:numbits + 54);
N = floor(length(x) / sstep);
encoded = LDPC_encoder(x,H,N);

interleaved = interleaver(encoded);
symbols_ak = bitmap(interleaved.').';

% save('Input_symbols.mat', 'symbols_ak', 'x', 'H', 'sstep');