clc; close all; clear global; clearvars;

% Initial bits
L = 2^20-1;
x = [PNSeq(L); PNSeq(L)];
% Matlab LDPC Encoder
H = comm.LDPCEncoder();         % parity check matrix

% encode block of 32400 bits
sstep = 32400;
% avoid final block to have less bit then required
numbits = floor(length(x) / sstep) * sstep;
x = x(1:numbits + 54);
% number of encoded packets
N = floor(length(x) / sstep);
% LDPC encoder
encoded = LDPC_encoder(x,H,N);
% interleaver
interleaved = interleaver(encoded);
% QPSK modulation with gray coding
symbols_ak = bitmap(interleaved.').';

save('Input_symbols.mat', 'symbols_ak', 'x', 'H', 'sstep');