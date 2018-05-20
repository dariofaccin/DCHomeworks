function Hd = GAA_filter
% All frequency values are normalized to 1.
Fpass = 0.45;           % Passband Frequency
Fstop = 0.55;           % Stopband Frequency
Dpass = 0.05;			% Passband Ripple
Dstop = 0.01;           % Stopband Attenuation
dens  = 20;             % Density Factor
[N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);
g_AA  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(g_AA);
save('GAA_filter.mat');