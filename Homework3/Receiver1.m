clc; close all; clear global; clearvars;

% Input signal
load('rec_input.mat');
T = 1;
L = 1023;
sigma_a = 2;

snr_db = 10;
snr_lin = 10^(snr_db/10);

% Matched filter
gm = qc(end:-1:1);

% h
h = conv(gm, qc);
% sample h with period T
h_T = downsample(h,4);
t0 = find(h==max(h));

figure()
subplot(121), stem([-t0+1:1:length(h)-t0],h,'b');
% xlabel('nT/4'), xlim([0 length(h)-1]), grid on
subplot(122), stem([-4:1:4],h_T,'Color','red');
% xlabel('nT'), xlim([0 length(h_T)-1]), grid on

t0 = find(h==max(h));
w_tilde = conv(wc,gm);

%% Adaptive LE pag. 639
N1 = 2;          % precursors 
N2 = 2;          % postcursors
M2 = 0;
M1 = 30;
packet = in_bits;
D = 5;
est_sigmaw = sigma_w;

% Zero padding of the i.r.
nb0 = 60;
nf0 = 60;
h_zero_ind = nb0+ N1 + 1;
hi = [zeros(nb0,1); h; zeros(nf0,1)];

p = zeros(M1, 1);
for i = 0:(M1 - 1)
    p(i+1) = sigma_a * conj(hi(N1+nb0+1+D-i));
end
R = zeros(M1);
for row = 0:(M1-1)
    for col = 0:(M1-1)
        first_sum = (hi((nb0+1):(N1+N2+nb0+1))).' * ...
            conj(hi((nb0+1-(row-col)):(N1+N2+nb0+1-(row-col))));
        second_sum = (hi((N1+nb0+1+1+D-col):(N1+nb0+1+M2+D-col))).' * ...
            conj((hi((N1+nb0+1+1+D-row):(N1+nb0+1+M2+D-row))));
        r_w = (row == col) * est_sigmaw; % This is a delta only if there is no g_M.
        
        R(row+1, col+1) = sigma_a * (first_sum - second_sum) + r_w;
        
    end
end

c_opt = inv(R)*p








% y = filter(q_mf,1,r_c);
% y = y(t0:end);
% y = downsample(y,4);
% detected = zeros(length(y),1);
% 
% for i=1:length(detected)
%     detected(i) = QPSK_detector(y(i));
% end
% numerrs = 0;
% for i=t0:length(detected)
%     if ( detected(i) ~= in_bits(i-t0+1) )
%         numerrs = numerrs + 1;
%     end
% end






