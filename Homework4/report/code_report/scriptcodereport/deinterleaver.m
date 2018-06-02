function [deinterleaved_bits] = deinterleaver(bits)
% This function receives a sequence of bits and unscrambles it
deinterleaved_bits = zeros(1,length(bits));

% The deinterleaver is just an interleaver with rows and cols switched
rows = 42;
columns = 42;

% We work with a rowsxcolumns matrix
for matrix = 0:(length(bits)/(rows*columns) - 1)
    curr_matrix = matrix * rows * columns;
    for col = 0:(columns-1)
        deinterleaved_bits(curr_matrix + col * rows + 1 : curr_matrix + col * rows + rows) = ...
            bits(curr_matrix + col + 1 : columns : curr_matrix + col + columns * rows);
    end
end
end