function mySystem = updateSystem(mySystem, TUcells, Fcells, M1cells, M2cells, TUprop,M1prop, M2prop, Fprop, L, Lt, i)

%Copy variables back to mySystem
mySystem.TU.TUcells = TUcells;
mySystem.TU.TUprop.isStem = TUprop.isStem;
mySystem.TU.TUprop.Pcap = TUprop.Pcap;
mySystem.TU.TUprop.isRes = TUprop.isRes;

mySystem.F.Fcells = Fcells;
mySystem.F.Fprop.Pcap = Fprop.Pcap;

mySystem.M1.M1cells = M1cells;
mySystem.M1.M1prop.Kcap = M1prop.Kcap;
mySystem.M1.M1prop.engaged = M1prop.engaged;

mySystem.M2.M2cells = M2cells;
mySystem.M2.M2prop.Kcap = M2prop.Kcap;
mySystem.M2.M2prop.engaged = M2prop.engaged;

mySystem.grid.L = L;
mySystem.grid.Lt = Lt;

if isfield(mySystem.grid, 'StepsDone')
    mySystem.grid.StepsDone = mySystem.grid.StepsDone + 1;
else 
    mySystem.grid.StepsDone = i;
end

end