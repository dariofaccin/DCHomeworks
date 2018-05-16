function [decisions] = equalization_pointC(x, c, b, D)
%EQUALIZATION for DFE
M2 = length(b);
y = conv(x,c);
y = downsample(y,2);
y = y(1:floor(length(x)/2));
detected = zeros(ceil(length(x)/2) + D, 1); 
for k=0:length(y)-1
     if (k <= M2)
        a_past = [flipud(detected(1:k)); zeros(M2 - k, 1)];
    else
        a_past = flipud(detected(k-M2+1:k));
    end
detected(k + 1) = QPSK_detector(y(k + 1) + b.'*a_past);
end
%scatterplot(y)
decisions = detected(D + 1:end);
end
