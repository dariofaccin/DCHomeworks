function [ Pe, count_err ] = SER( sent, detected )
% Computes the BER, it accepts symbols

if (length(sent) ~= length(detected))
    disp('Error in IBMAP, the sequences do not have the same length')
    return
end

count_err = 0;
for i=1:length(sent)
    if sent(i)~=detected(i)
        count_err = count_err + 1;
    end
end


Pe = count_err/length(sent);

end
