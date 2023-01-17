function [mySystem,cnst,valVars,fvals]=ABM_PSO_function(t,m,gridsize,nSteps,varNames,expname,lb,ub,TUcellNo,M1cellNo,M2cellNo,FcellNo,version,name,optno)

%Run particile swarm optimization 

addpath('./functions/') %Include all functions in path.
addpath('./data/')
addpath('./subroutines_ABM/')

opt = true;                                    %set if the optimization should run several times                                  

[mySystem,cnst] = getSystemParams(gridsize);   %get system parameters (168000 spots, start with 15000 cells, *11)
mySystem = getHyperParameters(expname,mySystem);
cnst.nSteps = nSteps;                               %set number of steps for the simulation (16*0.5(days)*24(hrs)=192 hours)
cnst.video = false;                             %don't show results on screen
cnst.simnumber = 1;
mySystem.params.TUcellNo = TUcellNo;            %set number of initial tumor cells
mySystem.params.M1cellNo = M1cellNo;            %number of initial M1 cells
mySystem.params.M1influxProb = 0;               %no macrophage influx
mySystem.params.M2cellNo = M2cellNo;            %number of initial M1 cells
mySystem.params.M2influxProb = 0;               %no macrophage influx
mySystem.params.FcellNo = FcellNo;                    %no fibroblasts
mySystem.params.TUps = 0;                       %no symmetric division as there are no stem cells in vitro

costFunction = @(x) MSE(x,mySystem, cnst,  varNames, version, t, m);   %create cost function for PSO

options = optimoptions('particleswarm','Display','iter');       %set some options for the PSO

%if the pso should be run several times, start loop, else, just run pso
%once
if opt       
    for i=1:optno
        [vals,fval,exitflag,output] = particleswarm(costFunction, length(varNames), lb, ub,options); %PSO
        fvals(i)=fval;
        valVars(i,:)=vals;
        disp("Finished simulation "+name+" "+num2str(i))
        save (name+"_valVars.mat",'valVars')
        save (name+"_fvals.mat",'fvals')
    end
else
    [valVars,fval,exitflag,output] = particleswarm(costFunction, length(varNames), lb, ub,options); %PSO
end

