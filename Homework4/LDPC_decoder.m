function [decoded] = LDPC_decoder(deinterleaved)

H = comm.LDPCDecoder('DecisionMethod','Hard decision');

numInfoBits = 32400;
decoded_bits = zeros(1,length(deinterleaved)/2);

% Iterate over the input info bits and decode them
for idx = 0:(length(deinterleaved)/(2*numInfoBits))-1
    current_bits = deinterleaved(2*idx*numInfoBits + 1 : 2*idx*numInfoBits + 2*numInfoBits);
    decoded_bits(idx*numInfoBits+1:idx*numInfoBits + numInfoBits) = step(H,current_bits.');
end
decoded = decoded_bits;
end