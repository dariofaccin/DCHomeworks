function [interleaved_bits] = interleaver(input)
    
    interleaved_bits = zeros(1,length(input));
    
    rows = 42;
    columns = 42;
    
    % We work with a rowsxcolumns matrix
    for matrix = 0:(length(input)/(rows*columns) - 1)
        curr_matrix = matrix * rows * columns;
        for col = 0:(columns-1)
            interleaved_bits(curr_matrix + col * rows + 1 : curr_matrix + col * rows + rows) = ...
                input(curr_matrix + col + 1 : columns : curr_matrix + col + columns * rows);
        end
    end
end