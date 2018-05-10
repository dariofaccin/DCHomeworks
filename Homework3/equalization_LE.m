function [decisions] = equalization_LE(x, c, M1, D, norm_fact)
%EQUALIZATION for LE 

y = zeros(length(x) + D , 1); % output of ff filter
detected = zeros(length(x) + D, 1); % output of td

for k = 0:length(x) - 1 + D
    
    if (k < M1 - 1)
        xconv = [flipud(x(1:k+1)); zeros(M1 - k - 1, 1)];
        
    elseif k > length(x)-1 && k < length(x) - 1 + M1
        xconv = [zeros(k-length(x)+1, 1); flipud(x(end - M1 + 1 + k - ...
            length(x) + 1:end) ) ];
        
    elseif k >= length(x) - 1 + M1 % just in case D is greater than M1
        xconv = zeros(M1, 1);
        
    else
        xconv = flipud(x(k - M1 + 1 + 1:k + 1));
    end 
    y(k+1) = c.'*xconv / norm_fact;
    detected(k + 1) = QPSK_detector(y(k + 1));
end
% scatterplot(y);
decisions = detected(D + 1:end);
end
