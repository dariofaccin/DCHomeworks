function [coded] = LDPC_encoder(input,h, N)

packets_num = N;
packet_length = 32400;
coded = step(h,input);
% for i=0:packets_num-1
%     uncoded = input(i*packet_length+1:i*packet_length+packet_length);
%     coded(:,i+1) = step(h,uncoded);
% end
% coded = reshape(coded,N*packet_length*2,1);

end