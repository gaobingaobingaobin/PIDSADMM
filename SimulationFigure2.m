
clc
clear
close all


%% TV data which extracted from simulationTable2.m
tempY{1}=[214,185;260,209;231,206;57,44;121,106];
tempY{2}=[176,154;213,174;189,171;47,37;99,88];
for i=1:2% 不同的tolerance 
subplot(1,2,i)

bar(tempY{i},1)

xlabel('Dimension of $D$','Interpreter','latex','fontsize',12);
ylabel('Iteration Numbers')
legend('PD-SADMM','PID-SADMM')
set(gca,'XTickLabel',{'100*100','200*200','300*300','400*400','500*500'})
end
