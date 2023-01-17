%%
%Grow castration resistant tumor with varying number of cell types.

%Clear variables
clear all;close all;clc
addpath('./functions/') %Include all functions in path.
addpath('./data/')
addpath('./subroutines_ABM/')
%%
%Define number of cells in simulations (McellNos(1) will be combined with
%FcellNos(1)).
McellNos = [125 500 500];
FcellNos = [125,125,500];

simnumber = 50; %How many simulations are run to generate median plot
j=1; %j = 1 M1 macrophage, j = 2 M2 macrophage, else = no macrophages
vid = false; %if vid is true, simulation will be shown and saved (advised to only set to true if simnumber = 1).

for i = 1:length(McellNos)
    McellNo = McellNos(i);
    FcellNo = FcellNos(i);
    
    %Define system
    [mySystem,cnst] = getSystemParams([125,125]);
    mySystem = getHyperParameters('Castration_Resistance',mySystem);
    cnst.nSteps = 60;
    cnst.video = vid;
    mySystem.params.TUpres = 0;
    mySystem.params.TUcellNo = 500;
    mySystem.params.TUresNo = 5;
    mySystem.params.FcellNo = FcellNo;
    
    %Select macrophage type & create name of simulation
    if j== 1
        mySystem.params.M1cellNo = McellNo;
        name = "CRPC-seeded-M1-"+num2str(McellNo/500)+"-1-"+num2str(FcellNo/500);
    elseif j == 2 
        mySystem.params.M2cellNo = McellNo;
        name = "CRPC-seeded-M2-"+num2str(McellNo/500)+"-1-"+num2str(FcellNo/500);
    else
        name = "CRPC-seeded-"+num2str(McellNo/500)+"-1-"+num2str(FcellNo/500);
    end
    
    %run simulation with current parameters for simnumber times
    for k = 1:simnumber 
        [mySystemn,TUcellNo,~,finalImage] = growTumor(mySystem,cnst);
        TUnoRelative(:,k,i) = (TUcellNo(:,1))/(TUcellNo(1,1));   
        TUresRelative(:,k,i) = (TUcellNo(:,2))/(TUcellNo(1,2));
    end
    
    %logtransform cell numbers 
    logTUnoRelative = log(TUnoRelative(:,:,i)+1);
    logTUresRelative = log(TUresRelative(:,:,i)+1);
    TUnoRelativemed = median(logTUnoRelative,2);
    
    %calculate iqrange of simulation outputs & set up some parameters for
    %plotting
    TUnoRelativeq = quantile(logTUnoRelative,[0.25 0.75],2);
    TUnoRelativeq = [TUnoRelativeq(:,2)-TUnoRelativemed,TUnoRelativemed-TUnoRelativeq(:,1)];
    TUresRelativemed = median(logTUresRelative,2);

    TUresRelativeq = quantile(logTUresRelative,[0.25 0.75],2);
    TUresRelativeq = [TUresRelativeq(:,2)-TUresRelativemed,TUresRelativemed-TUresRelativeq(:,1)];
    
    %Create plots of median & iqrange TUnoRelative over time 
    figure()
    time = 0:4:cnst.nSteps*4;    
    shadedErrorBar(time,TUnoRelativemed,TUnoRelativeq,'lineProps',{'color', 'red', 'linewidth',2})
    hold on
    shadedErrorBar(time,TUresRelativemed,TUresRelativeq,'lineProps',{'color', 'black', 'linewidth',2})
    xlabel('Time (hrs)'); 
    ylabel('log(Relative tumor cell number+1)')
    title(name)
    set(gca, 'FontSize', 15)
    ylim([0,6])
    legend(["Normal Tumor Cells","Resistant Tumor Cells"],'Location','NorthWest')
    
    %save video or save plots
    if vid == true
         writeMyVideo(finalImage,"output/video-"+name,4)
    else
        saveas(gcf,"output/plot-median-"+name+"-median.pdf")
    end
    
end
%%
