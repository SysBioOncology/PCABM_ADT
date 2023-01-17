function plotPSOVar3D(vars,fvals,names,labels)
%Scatter plot optimized variables 2D

figure()

for i=1:size(fvals,1)
    fval=fvals(i,:);
    valVars=vars(:,:,i);
    numtrail=length(fval);
    
    scatter3(valVars(:,1),valVars(:,2),valVars(:,3),'filled')
    hold on
end

xlim([0,Inf])
xlabel(labels(1))
ylim([0,Inf])
ylabel(labels(2))
zlim ([0,Inf])
zlabel(labels(3))
title('Distribution of optimized parameters')
legend(names)
set(gca,'FontSize',15)

end
