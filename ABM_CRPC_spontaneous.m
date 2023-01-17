%Simulate spontaneous development of CRPC over time.
%Starts with normal tumor cell growth with R1881, system gets castrated and
%switches to DMSO conditions. After a certain amount of time, tumor cells
%can become resistant upon proliferation. 

%Clear variables
clear all;close all;clc
%Include all necessary directories in path
addpath('./functions/') 
addpath('./data/')
addpath('./subroutines_ABM/')

vidname = "output/CRPC_full"; %Name of simulation

%Get system template
[mySystem,cnst] = getSystemParams([125,125]);   
%Choose number of cells in the system
mySystem.params.TUcellNo = 500;
mySystem.params.M1cellNo = 0;
mySystem.params.M2cellNo = 0;
mySystem.params.FcellNo = 0;
cnst.video = true;

%Start with growing cells 'normally' in presence of hormone
mySystem1 = getHyperParameters('R1881',mySystem);
cnst1 = cnst;
cnst1.nSteps = 40;

[mySystem2,~,~,finalImage1] = growTumor(mySystem1,cnst1);

%Use template of normally grown cells and remove hormone
%Change to DMSO condition
mySystem2 = getHyperParameters('DMSO',mySystem2);
cnst2 = cnst; 
cnst2.nSteps = 40;
cnst2.newSystem = false;

[mySystem3,~,~,finalImage2] = growTumor(mySystem2,cnst2);

%Use template of cells grown without hormone
%Move to Castration Resistance
mySystem3 = getHyperParameters('Castration_Resistance',mySystem3);
cnst3 = cnst;
cnst3.nSteps = 100;
cnst3.newSystem = false;

[mySystem4,~,~,finalImage3] = growTumor(mySystem3,cnst3);

%Write Video
finalImage =[finalImage1,finalImage2,finalImage3];
writeMyVideo(finalImage,vidname,4)

