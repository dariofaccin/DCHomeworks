function [coded] = LDPC_encoder(input, h, N)

packets_num = N;
packet_length = 32400;
coded = zeros(length(input)*2,1);
for i=0:packets_num-1
    uncoded = input(i*packet_length+1:i*packet_length+packet_length);
    coded(2 * i * packet_length + 1:2 * i * packet_length + 2 * packet_length) = step(h, uncoded);
end

coded = reshape(coded,length(input)*2,1);

end