
clc
clear
close all
%% LASSO data which extracted from simulationTable1.m
tempY{1}=[35,30;36,32;33,28;36,30;30,26];
tempY{2}=[34,31;35,32;32,29;35,31;28,25];

for i=1:2% 

    
subplot(1,2,i)

bar(tempY{i},1)

xlabel('Dimension 0f $A$','Interpreter','latex','fontsize',12);
ylabel('Iteration Numbers')
legend('PD-SADMM','PID-SADMM')
set(gca,'XTickLabel',{'900*3000','1050*3500','1200*4000','1350*4500','1500*5000'})

end
