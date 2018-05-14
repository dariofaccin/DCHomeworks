clc; close all; clear global; clearvars;

load('Viterbi_out.mat');

Pe = zeros(8,1);

printmsg_delete = '';

for i=1:8
    printmsg = sprintf('iteration = %d\n', i);
    fprintf([printmsg_delete, printmsg]);
    printmsg_delete = repmat(sprintf('\b'), 1, length(printmsg));
    nerr = length(find(in_bits(i:i+length(detected)-1)~=detected));
    Pe(i,1) = nerr/length(in_bits(i:i+length(detected)-1));
end