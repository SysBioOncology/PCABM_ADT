function plotPSOVarCat(vars,labels,varname)
%Create a categorical scatter plot of optimized variables (1D)

figure()

for i=1:size(vars,2)
    var=vars(:,i);
    xdata=repmat(i,length(var),1);
    scatter(xdata(:), var, 'filled', 'jitter','on', 'jitterAmount', 0.25);
    hold on
end

xlim([0,size(vars,2)+1])
xticks(0:1:size(vars,2)+1)
xticklabels([{' '},labels,{' '}])
%ylim([0,0.01])
ylim([0,Inf])
ylabel(varname)
title('Distribution of optimized parameters')
set(gca,'FontSize',15)

end
