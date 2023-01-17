function mse = MSE(x,mySystem, cnst,  varNames, version, t, m)
%Calculate MSE for PSO

time = 0:4:cnst.nSteps*4;
mySystemO = mySystem;

%Adjust mySystem parameters for input parameters in x
for i = 1:length(varNames)
    varName = varNames(i);
    val = x(i);
end

for j = 1:cnst.simnumber
    [~,TUcellNo,FcellNo,~] = growTumor(mySystem,cnst);                %run ABM
    TUcellNo=TUcellNo(:,1);
    if version == "TU"
        CellnoRelative = 1+((TUcellNo-TUcellNo(1))/(TUcellNo(1)));    %calculate relative cell number either tumor cells "TU" or fibroblasts "F"
    elseif version == "F"
        CellnoRelative = 1+((FcellNo-FcellNo(1))/(FcellNo(1)));
    else
        disp("Wrong version selected")
    end
    CellnoInt (:,j) = interp1(time,CellnoRelative,t);               %interpolate relative cell number to match time steps of incucyte
end

CellnoMean= mean(CellnoInt,2);                                   %calculate mean relative cell number


mse = mean((m - CellnoMean).^2);                               %calculate MSE between incucyte and model output

end