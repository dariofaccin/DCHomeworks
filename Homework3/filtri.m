%% C
close all
figure, stem(0:length(c_opt)-1,abs(c_opt)), hold on, grid on
ylabel('$|c|$'), xlabel('$n\frac{T}{2}$'); xlim([0 length(c_opt)-1]);
%% PSI
close all
% figure, stem(-2:length(psi)-3,abs(psi)),
stem(-(find(psi==max(psi))-D)+1:length(psi)-(find(psi==max(psi))-D+1)+1,abs(psi))
xlim([-(find(psi==max(psi))-D)+1 length(psi)-(find(psi==max(psi))-D+1)+1]), grid on
ylabel('$|\psi|$'), xlabel('$n\frac{T}{2}$'); 
%% B
close all
figure, stem(0:length(b)-1,abs(b)), hold on, grid on
ylabel('$|b|$'), xlabel('n'); xlim([0 length(b)-1]);
%% gm
close all
figure, stem(abs(gm));
xlabel('$m\frac{T}{4}$');
ylabel('$g_m$')
xlim([1 length(gm)]);
grid on


stem(-(find(psi==max(psi))-D)+1:length(psi)-(find(psi==max(psi))-D+1)+1,abs(psi))