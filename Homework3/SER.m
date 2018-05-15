function [Pe, count_errors] = SER(sent, detected)
% Computes the symbol-error rate, it accepts QPSK symbols
count_errors = 0;
for i=1:length(sent)
    if sent(i) ~= detected(i)
        count_errors = count_errors + 1;
    end
end
% count_errors = sum((sent-detected)~=0);
Pe = count_errors/length(sent);
end