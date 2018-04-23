function [welch_est, Ns] = welchPSD(inputsig, window, overlaps)

D = length(window);             % Length of the window
K = length(inputsig);           % Length of input signal
Mw = sum(window .^ 2) * (1/D);  % Normalized energy of the window
N_s = floor((K-D)/(D-overlaps) + 1);        % Number of subsequences
P_per = zeros(K, N_s);          %Initialization of each periodogram

for s = 0:(N_s-1)
    x_s = window .* inputsig(s*(D-overlaps)+1:s*(D-overlaps)+D);
    X_s = fft(x_s, K);
    P_per(:,s+1) = (abs(X_s)).^2 * (1/(D*Mw)); % Periodogram for the window
end
welch_est = sum(P_per, 2) * (1/N_s);        % Sum of all periodograms
Ns = length(welch_est);
end