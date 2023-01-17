function [mySystem, cnst]=getSystemParams(space)
%Sets the parameters needed for modeling the ABM (standard is condition
%with R1881)
mySystem.params.initialSeed = 1;   % initial random seed, default 1

%DEFINE TUMOR CELL PARAMETERS
mySystem.params.TUpprol = 0.0846;  %probability of proliferation
mySystem.params.TUpmig = 0.1167;     %probability of migration
mySystem.params.TUpdeath = 0.00284;%probability of death
mySystem.params.TUrwalk = 0.5;    %random influence on movement
mySystem.params.TUpmax = 4;   %initial proliferation capacity
mySystem.params.TUps = 0;          %probability of symmetric division
mySystem.params.TUcellNo = 1500;      %center first tumor cells or use multiple randomly distributed
mySystem.params.TUresNo = 0;        %seed some resitant cells
mySystem.params.TUpres = 0;        %probability of proliferation leading to resistant tumor cell 
mySystem.params.TUpprolres = 0;    %probability of proliferation of resistant cells
mySystem.params.TUpmigres = 0;     %migration of resistant cells
mySystem.params.TUpmaxres = 0;      %poliferation capacity of resistant cells

%DEFINE MACROPHAGE TYPE 1 PARAMETERS
mySystem.params.M1kmax = 11;        %killing capacity 
mySystem.params.M1pkill = 0.0306;   %probability of killing 
mySystem.params.M1pmig = 0.2667;    %probability of migration
mySystem.params.M1pdeath = 0.0049;  %probability of death
mySystem.params.M1rwalk = 0.8;      %random influence on movement
mySystem.params.M1speed = 40;       %speed of movement
mySystem.params.M1influxProb = 0; %probability of influx
mySystem.params.M1influxRate = 1;   %number of cell influx
mySystem.params.M1cellNo = 0;       %number of initial cells
mySystem.params.M1engagementDuration = 60; %number of steps immune cell is engaged


%DEFINE MACROPHAGE TYPE 2 PARAMETERS
mySystem.params.M2kmax = 11;   %killing capacity 
mySystem.params.M2pkill = 0.0127;   %probability of killing 
mySystem.params.M2pkill = 0;        %probability of cell promotion
mySystem.params.M2pmig = 0.2667;    %probability of migration
mySystem.params.M2pdeath = 0.0049;  %probability of death
mySystem.params.M2rwalk = 0.8;      %random influence on movement
mySystem.params.M2speed = 40;       %speed of movement
mySystem.params.M2influxProb = 0; %probability of influx
mySystem.params.M2influxRate = 1;   %number of macrophages each influx
mySystem.params.M2cellNo = 0;       %number of initial cells
mySystem.params.M2TUadd = 0;         %Tumor promoting
mySystem.params.M2engagementDuration = 60; %number of steps immune cell is engaged

%DEFINE FIBROBLAST PARAMETERS
mySystem.params.Fpprol = 0.0838;      %probability of proliferation
mySystem.params.Fpmig = 0.4;       %probability of migration
mySystem.params.Fpdeath = 0.0018;    %probability of death
mySystem.params.Fpmax = 4;          %initial proliferation capacity
mySystem.params.Frwalk = 0.5;      %random influence on movement
mySystem.params.FcellNo = 0;     %ratio of initial number of cells

%DEFINE DIMENSION PARAMETERS
%to create NxM grid
mySystem.grid.N=space(1);
mySystem.grid.M=space(2);

%DEFINE SOME CONSTANTS
cnst.nSteps = 35;                    %number of steps in the simulation
cnst.video = false;                    %draw video is true
cnst.newSystem = true;              %initialize new system or use previously defined
