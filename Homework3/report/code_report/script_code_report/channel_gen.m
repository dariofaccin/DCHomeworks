clc; close all; clear global; clearvars;
set(0,'defaultTextInterpreter','latex')

load('Useful.mat');
T = 1;              % Symbol period
Tc = T/4;           % upsampling period
Q = T/Tc;           % Interpolation factor
snr_db = 10;
snr_lin = 10^(snr_db/10);
sigma_a = 2;		% Input variance

L = 1023;			% Input signal: from a PN sequence generate QPSK
x = PNSeq(L);
in_bits = bitmap(x(1:length(x)-1));