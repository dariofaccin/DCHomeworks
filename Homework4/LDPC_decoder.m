function [decoded] = LDPC_decoder(deinterleaved)

H = comm.LDPCDecoder('DecisionMethod','Hard decision');

packets_num = 32;
packet_length = 64800;
% decoded = zeros(packet_length/2,packets_num);
% for i=0:packets_num-1
%     temp = deinterleaved(i*packet_length+1:i*packet_length+packet_length);
%     decoded(:,i+1) = step(h,temp);
% end
% decoded = reshape(decoded,packets_num*packet_length/2,1);

numInfoBits = 32400;
decoded_bits = zeros(1,length(deinterleaved)/2);
% Iterate over the input info bits and decode them
% for idx = 0:(length(deinterleaved)/(2*numInfoBits))-1
%     current_bits = deinterleaved(2*idx*numInfoBits + 1 : 2*idx*numInfoBits + 2*numInfoBits);
%     decoded_bits(idx*numInfoBits+1:idx*numInfoBits + numInfoBits) = step(H,current_bits);
% end
decoded_bits = step(H,deinterleaved);
decoded = decoded_bits;


end




