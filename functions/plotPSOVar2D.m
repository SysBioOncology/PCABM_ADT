function plotPSOVar(vars,names,labels)
%Scatter plot optimized variables 2D

figure()

for i=1:2:(size(vars,2)-1)
    scatter(vars(:,i),vars(:,i+1),'filled')
    hold on
end

hold on
xlim([0,Inf])
xlabel(labels(1))
ylim([0,Inf])
ylabel(labels(2))
title('Distribution of optimized parameters')
legend(names,'Location','NorthWest')
set(gca,'FontSize',15)

end

