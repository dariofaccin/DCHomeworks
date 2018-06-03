%% square root rised cosine centered
stem(-(length(g_srrc)-1)/2:(length(g_srrc)-1)/2,g_srrc)
xlim([-(length(g_srrc)-1)/2 (length(g_srrc)-1)/2 ]), grid on
xlabel('$nT_c$')
ylabel('$g_{\sqrt{rcos}}(nT_c)$')

%% square root rised cosine not centered
stem(g_srrc)
xlim([1 length(g_srrc)-1]), grid on
xlabel('$nT_c$')
ylabel('$g_{\sqrt{rcos}}(nT_c)$')

%% square-root raised cosine in f
[G_srrc, f] = freqz(g_srrc,1,1024,'whole');
f = f/(2*pi);
figure, plot(f,10*log10(abs(G_srrc))), grid on, xlim([0 0.5])
ylim([-15 5])
xlabel('$f/T_c$')
ylabel('$|G{\sqrt{rcos}}|$ [dB]');

%% qc from homework 3
figure, stem(0:length(qc)-1,qc)
xlabel('$nT_c$')
ylabel('$q_c$')
xlim([0 length(qc)-1]), grid on

%% qr
figure, stem(q_r), grid on
xlabel('$nT_c$')
ylabel('$q_r(nT_c)$')
xlim([1 length(q_r)-1]);

%% Qr
[Q_r, f] = freqz(q_r,1,1024,'whole');
f = f/(2*pi);
plot(f,10*log10(abs(Q_r)))
grid on
xlim([0 0.5]);
ylim([-5 10]);
xlabel('$f/T_c$')
ylabel('$|Q_R|$ [dB]');

%% h(mT_OFDM)
figure
peak_pos = find(q_r == max(q_r));
stem(-(peak_pos-1):length(q_r)-(peak_pos-1)-1,q_r)
q_r_ds = q_r(1:4:end);
stem(-(peak_pos-1)/4:length(q_r_ds)-(peak_pos-1)/4-1,q_r_ds)
xlim([-(peak_pos-1)/4 length(q_r_ds)-(peak_pos-1)/4-1])
grid on
xlabel('$mT_{OFDM}$');
ylabel('h')

%% DFT of h
H = fft(q_r_ds,512);
f = linspace(0,1,512);
plot(f,10*log10(abs(H/length(H))))
xlim([0 1]);
xlabel('f$\mathcal{M}/T_{OFDM}$')
ylabel('H [dB]')
grid on

