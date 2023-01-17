function [CellnoRelativeFinal,miniq,maxiq,vars]=iqsensanalysis(vars,varnames,TUcellNoinput,M1cellNo,M2cellNo,FcellNoinput,version,name,hyperparams)
%Run sensitivity analysis for sets of optimized variables

%Set system parameters
[mySystem,cnst] = getSystemParams([125,125]);
mySystem = getHyperParameters(hyperparams,mySystem);
mySystem.params.TUcellNo=TUcellNoinput; 
mySystem.params.M1cellNo=M1cellNo;
mySystem.params.M2cellNo=M2cellNo;
mySystem.params.FcellNo=FcellNoinput;

time=0:4:cnst.nSteps*4;

%Grow tumor with median parameters
[~,TUcellNo,FcellNo,~] = growTumor(mySystem,cnst);
if version =="TU"
    CellnoRelativeFinal = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); %calculate relative cell number 
else
    CellnoRelativeFinal = 1+((FcellNo-FcellNo(1))/(FcellNo(1)));
end

q=quantile(vars,[0.25 0.75]); %Get interuquantile range
sizeq=size(q);

%Calculate tumor cell growth for parameters in interquantile range 
for j =1:sizeq(1)
    
    [mySystem,cnst] = getSystemParams([125,125]);
    mySystem = getHyperParameters(hyperparams,mySystem);
    mySystem.params.TUcellNo=TUcellNoinput;
    mySystem.params.M1cellNo=M1cellNo;
    mySystem.params.M2cellNo=M2cellNo;
    mySystem.params.FcellNo=FcellNoinput;

     varname=varnames(j);
    
    if sizeq(1) == 1
        varmin=q(1);
        varmax=q(2);
    else
        varmin=q(1,j);
        varmax=q(2,j);
    end
      
    %Remove variables outside of iq range
    for i=1:length(vars)
        if (vars(i,j)<varmin) || (vars(i,j)>varmax)
        vars(i,j)=nan;
        end
    end
   
    %Create range of 20 variables to vary per parameter
    varssamples = linspace(varmin,varmax,20);
    
    for i=1:length(varssamples)
        mySystem.params.(varname)=varssamples(i);
        [~,~,TUcellNo, ~, FcellNo,  ~,~] = growTumor(mySystem,cnst);
        
        if version == "TU"
            CellnoRelativeVar(:,i)=1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); 
        else
            CellnoRelativeVar(:,i)=1+((FcellNo-FcellNo(1))/(FcellNo(1)));
        end
    end
 
end

%Get set of random variables from actual optimized variables within
%quantile

varssamples=getsets(vars);

for i = 1:length(varssamples)
    varssample=varssamples(i,:);
    
    for j=1:length(varssample)
        mySystem.params.(varnames(j))=varssample(j);
    end
    
    [~,~,TUcellNo, ~, FcellNo,  ~,~] = growTumor(mySystem,cnst);
    
    if version == "TU"
        CellnoRelativeVars(:,i)=1+((TUcellNo-TUcellNo(1))/(TUcellNo(1))); 
    else
        CellnoRelativeVars(:,i)=1+((FcellNo-FcellNo(1))/(FcellNo(1)));
    end
end

    %Range of simulations in iqrange from median
    [m,i]=max(CellnoRelativeVars(end,:));
    maxiq=CellnoRelativeVars(:,i);
    for i=1:length(maxiq)
        if CellnoRelativeFinal(i) > maxiq(i)
            maxiq(i) = CellnoRelativeFinal(i);
        end
    end
    maxiq=maxiq'-CellnoRelativeFinal;
    [m,i]=min(CellnoRelativeVars(end,:));
    miniq=CellnoRelativeVars(:,i);
    miniq=CellnoRelativeFinal-miniq';

end
