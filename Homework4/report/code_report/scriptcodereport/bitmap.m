function [output] = bitmap(input)
% Check if the input array has even length
L = length(input);

output = zeros(1,L);

% Map each couple of values to the corresponding symbol
for idx = 1:2:L-1
    if (isequal(input(idx:idx+1), [0; 0] ))
        output(idx) = -1-1i;
    elseif (isequal(input(idx:idx+1), [1; 0] ))
        output(idx) = 1-1i;
    elseif (isequal(input(idx:idx+1), [0; 1] ))
        output(idx) = -1+1i;
    elseif (isequal(input(idx:idx+1), [1; 1] ))
        output(idx) = +1+1i;
    end
end

output = output(1:2:end);
end