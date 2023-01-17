function plot_incucyte(filename,sheet,tit,names,ylab,colors, marker) 

%plot incucyte data

file = readtable(filename,'Sheet',sheet,'PreserveVariableNames',true);         %load dataset

figure()                            %create new figure

i=1;
c=1;
while i < width(file)
    
    co=colors(c);
    ma=marker(c);
    
    tinytable = file(:,i:i+2);      %get t,m and sd for one experiment
    x= tinytable{:,1};
    y= tinytable{:,2};
    error= tinytable{:,3};
    
    name = string(tinytable.Properties.VariableNames(1)); %get the name of the experiment
    
    errorbar(x,y,error,ma,'DisplayName',name,'Linewidth',1,'MarkerFaceColor',co,'MarkerEdgeColor',co,'Color',co) %plot the experiment
    hold on 
    
    i=i+3;
    c=c+1;
end

if isempty(names)
    legend
else
    legend(names)
end

xlabel("Time (hrs)")
ylabel(ylab)
title(tit)
set(gca,'FontSize',15)
ylim([0,9])

end