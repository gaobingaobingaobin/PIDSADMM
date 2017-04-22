
clc
clear
close all
%% LASSO data which extracted from simulationTable1.m
tempY{1}=[35,30;36,32;33,28;36,30;30,26];
tempY{2}=[34,31;35,32;32,29;35,31;28,25];

for i=1:2%
    
    
%     subplot(1,2,i)
    figure(i);
    b=bar(tempY{i},1);
    hold on;
    color=[0 0 0.75;0 1 0];
    for j=1:2
        set(b(j),'FaceColor',color(j,:));
    end
    
    xlabel('Dimension 0f $A$','Interpreter','latex','fontsize',12);
    ylabel('Iteration Numbers')
    legend('P-SADMM','ID-SADMM')
    set(gca,'XTickLabel',{'900*3000','1050*3500','1200*4000','1350*4500','1500*5000'})
    if i==1
        title('$\alpha=-0.3$','Interpreter','latex','fontsize',12);
    else
        title('$\alpha=0.3$','Interpreter','latex','fontsize',12);
    end
end
