function [L,Mcells,Mprop,TUcells,TUprop,Lt] = mCellRound(L,Lt,Mcells,TUcells, Mprop,TUprop,MinfluxProb,MinfluxRate,Mkmax,Mspeed,Mpmig,Mrwalk,Mpkill,Mpdeath,MengagementDuration,ChtaxMap,nh,cell3,cell4)

if MinfluxRate > 0            %Check if there even is an influxrate
if rand () <= MinfluxProb     %Check if there is an influx of cells
if sum(~L(:)) > 0             %Check if there are empty locations for the new cells
    [~,coordsNewMcells] = datasample(L(:),MinfluxRate,'Replace',false, 'Weights', uint8(~L(:))); %place new cells in empty location
    L(coordsNewMcells) = true;                              %place new cells on grid
    nNewCells = numel(coordsNewMcells);                     %number of new cells
    Mcells = [Mcells, coordsNewMcells];                     %add new cells to stack
    Mprop.Kcap = [Mprop.Kcap, repmat(Mkmax,1,nNewCells)];   %add properties
    Mprop.engaged = [Mprop.engaged, zeros(1,nNewCells)];    %add properties    
end
end
end

[Mcells,Mprop] = shuffleCells(Mcells,Mprop);    %shuffle cells
if numel(Mcells) > 0                            %if there are any macrophages
for j = 1:(Mspeed-1)                            %allow macrophages to move Mspeed-1 times per round
    L([Mcells,TUcells,cell3,cell4]) = true;                      %ensure new macrophages are on grid
    [L, Mcells] = M_go(Mcells, Mpmig, Mrwalk, ChtaxMap, L, nh);  %macrophage movement
    [TUcells, TUprop, Mcells, Mprop, L, Lt]  = M_kill (TUcells, TUprop, Mcells, Mprop, L, Lt, Mpkill, nh, ChtaxMap, MengagementDuration); %macrophage killing
    Mprop.engaged(Mprop.engaged > 0) = Mprop.engaged(Mprop.engaged > 0)-1; %un-engage macrophage
end

%allow macrophages to move once more or die
[L, Mcells, Mprop] = M_go_die(Mcells, Mprop, Mpmig, Mpdeath, Mrwalk, ChtaxMap, L, nh);


end
    
end