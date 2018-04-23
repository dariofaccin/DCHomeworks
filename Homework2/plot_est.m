function plot_est(cor, ls, sigdB)

set(0,'defaultTextInterpreter','latex')
% figure()
scrsz = get(0,'ScreenSize');
figure('Position',[15 scrsz(4)/5 scrsz(3)/1.5 scrsz(4)/1.5])
[N,lengL] = size(cor);
N = [1:1:N];
a = sigdB*ones(1,20); 

plot(N,ls(:,1),'b')
hold on, plot(N,ls(:,2),'c')
hold on, plot(N,ls(:,3),'g')
hold on, plot(N,ls(:,4),'y')
hold on, plot(N,ls(:,5),'m')
hold on, plot(N,ls(:,6),'r')
hold on, plot(N,cor(:,1),'bo--')
hold on, plot(N,cor(:,2),'co--')
hold on, plot(N,cor(:,3),'go--')
hold on, plot(N,cor(:,4),'yo--')
hold on, plot(N,cor(:,5),'mo--')
hold on, plot(N,cor(:,6),'ro--')

hold on, plot(a,'b--','LineWidth',2)

xlabel('$N_h$');
ylabel('$\mathcal{E}/L$ [dB]')
xlim([1 20]);
legend('3','L7','L15','L31','L63','L127') 
% legend('L7','L15','L31','L63','L127','255') 
set(gca,'FontSize',15);
grid on

end