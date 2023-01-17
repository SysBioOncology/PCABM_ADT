%Example of how to run the ABM

clear all;close all;clc
addpath('./functions/') %Include all functions in path.
addpath('./data/')
addpath('./subroutines_ABM/')

simname = "RUN"; %define name of video for saving

[mySystem,cnst] = getSystemParams([125,125]);   %get system parameters
mySystem = getHyperParameters('R1881',mySystem); %define condition
cnst.video = true;                             %show video on screen

%Adjust initial cell numbers (optional)
mySystem.params.TUcellNo = 1500;                
mySystem.params.M1cellNo = 50;
mySystem.params.M2cellNo = 50;
mySystem.params.FcellNo = 500;

%run model 
[mySystemnew,TUcells,~,finalImage] = growTumor(mySystem,cnst);

%Plot relative TUcellNo over time
TUcellNo = TUcells(:,1);
TUnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); %calculate relavite cell number 
time = linspace(0,cnst.nSteps*4,cnst.nSteps+1);

figure()
plot(time,TUnoRelative,'Color','k','LineWidth',3,'DisplayName','PCABM')
xlabel('Time (hrs)'); 
ylabel('Relative Number of Tumor Cells')
title("Simulated Tumor Cells over Time")
legend('Location','NorthWest')
set(gca, 'FontSize', 15)

writeMyVideo(finalImage,"output/vid_"+simname,4)
saveas(gcf,"output/plot_"+simname+".png")
