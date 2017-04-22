
clc
clear
close all


%% TV data which extracted from simulationTable2.m
tempY{1}=[214,185;260,209;231,206;57,44;121,106];
tempY{2}=[176,154;213,174;189,171;47,37;99,88];
for i=1:2% 涓嶅悓鐨則olerance 
% subplot(1,2,i)
 figure(i);
b=bar(tempY{i},1);
hold on;
    color=[0 0 0.75;0 1 0];
    for j=1:2
        set(b(j),'FaceColor',color(j,:));
    end
xlabel('Dimension of $D$','Interpreter','latex','fontsize',12);
ylabel('Iteration Numbers')
legend('P-SADMM','ID-SADMM')
  if i==1
        title('$\alpha=-0.1$','Interpreter','latex','fontsize',12);
    else
        title('$\alpha=0.1$','Interpreter','latex','fontsize',12);
    end
set(gca,'XTickLabel',{'100*100','200*200','300*300','400*400','500*500'})
end
