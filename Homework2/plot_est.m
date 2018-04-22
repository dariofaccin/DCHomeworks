function plot_est(cor, ls, sigdB)

set(0,'defaultTextInterpreter','latex') 
figure()
[N,lengL] = size(cor);
N = [1:1:N];
a = sigdB*ones(1,20); [3 7 15 31 63 127]

plot(N,cor(:,1),'b--')
hold on, plot(N,cor(:,2),'c--')
hold on, plot(N,cor(:,3),'g--')
hold on, plot(N,cor(:,4),'y--')
hold on, plot(N,cor(:,5),'m--')
hold on, plot(N,cor(:,6),'r--')
hold on, plot(N,ls(:,1),'b')
hold on, plot(N,ls(:,2),'c')
hold on, plot(N,ls(:,3),'g')
hold on, plot(N,ls(:,4),'y')
hold on, plot(N,ls(:,5),'m')
hold on, plot(N,ls(:,6),'r')
hold on, plot(a,'b--','LineWidth',2)

xlabel('$N_h$');
xlim([1 20]);
legend('L3','L7','L15','L31','L63','L127')          

end