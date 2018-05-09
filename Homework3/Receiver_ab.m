clc; close all; clear global; clearvars;

% Input signal
load('rec_input.mat');
T = 1;
L = 1023;
sigma_a = 2;

snr_db = 10;
snr_lin = 10^(snr_db/10);

% Matched filter
gm = conj(qc(end:-1:1));

% h
h = conv(qc,gm);
h = h(h>max(h)/100);
h = h(3:end-2);

% sample h with period T
h_T = downsample(h,4);

% figure()
% subplot(121), stem(-8:8,h,'b');
% xlabel('nT/4'), xlim([-9 9]), grid on
% subplot(122), stem(-2:2,h_T,'Color','red');
% xlabel('nT'), xlim([-3 3]), grid on

r_c_prime = filter(gm,1,r_c);

t0_bar = find(h_T==max(h_T));             % timing phase
r_c_prime = r_c_prime(t0_bar:end);
x = downsample(r_c_prime,4);

r_gm = xcorr(gm,gm);                     % crosscorr ??
rw_tilde = sigma_w*r_gm;

% save('Receiver_ab.mat')

%% RECEIVER a
clc; close all; clear global; clearvars;
load('rec_input.mat')
load('Receiver_ab.mat')

M1 = 5;
M2 = 0;
D = 5;
[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

y = filter(c_opt,1,x);
detected = zeros(length(y),1);

scatterplot(y)

for i=1:length(detected)
    detected(i) = QPSK_detector(y(i));
end

detected = detected(D+1:end);

numerrs = 0;
for i=1:length(detected)
    if ( detected(i) == in_bits(i))
        numerrs = numerrs +1;
    end
end
numerrs

%% RECEIVER b
M1 = 5;
M2 = 3;
D = 2;
[c_opt, Jmin] = Adaptive_DFE(h_T, rw_tilde, sigma_a, M1, M2, D);

