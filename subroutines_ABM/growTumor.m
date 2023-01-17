function [mySystem,TUcellNo,FcellNo,finalImage] = growTumor(mySystem,cnst)
%Performs the agent based modeling.
%Inputs are two structures; the parameters and some global constants
%defined in getSystemParams.
%Ouptut is the system at final iteration, the number of tumor cells over time (two
%columns, non-resistant and resistant respectively), the number of fibroblasts over time and the final video. 

%BUILD ENVIRONMENT
%Throw parameters to workspace:
cellfun(@(f) evalin('caller',[f ' = mySystem.params.' f ';']), fieldnames(mySystem.params));
cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));

%INITIALIZE SYSTEM
if cnst.newSystem
    %Create empty system:
    L=false(N,M);                                                      %initialize grid
    L(1,:) = true; L(:,1) = true; L(end,:) = true; L(:,end) = true;     %set grid boundary to occupied

    %Add first tumor cell
    if TUcellNo == 1                            %center location of first TUcell if there is only 1 cell
        TUcells = int32(N*round(M/2)-round(N/2)); 
        TUprop.isStem=true;                     %set property of first tumor cell to stemcell
        TUprop.isRes=false;                     %set first cell to be non-resistant
        TUprop.Pcap=uint8(TUpmax);              %set property of first tumor cell to max proliferation
    else                                        %randomly distribute tumor cells if there are more than 1
        C = (find(L==false))';
        TUcells = int32(C(randperm(numel(C),TUcellNo+TUresNo))); 
        TUprop.isStem = true(1,TUcellNo+TUresNo); 
        TUprop.isRes = [true(1,TUresNo),false(1,TUcellNo)]; 
        TUprop.Pcap = [uint8(TUpmaxres*ones(1,TUresNo)),uint8(TUpmax*ones(1,TUcellNo))]; 
        TUcellNo(1,2) = TUresNo;
    end
    
    %Randomly add fibroblasts to grid depending on the fibroblast ratio

    C = (find(L==false))';                                           %check which spots are not taken
    Fcells = int32(C(randperm(numel(C),FcellNo)));                  %randomly allocate positions for the fibroblasts, not taking tumor position
    Fprop.Pcap = uint8(repmat(Fpmax,1,FcellNo));                    %set property of fibroblasts to max proliferation
       

    %Randomly add defined number of macrophages type 1 and 2 to grid
    C = (find(L==false))' ;                            %check again free spots
    M1cells= int32(C(randperm(numel(C),M1cellNo)));    %randomly allocate positions for M1s
    M1prop.Kcap = uint8(repmat(M1kmax,1,M1cellNo));    %add properties: max killing capacity
    M1prop.engaged = uint8(zeros(1,M1cellNo));         %add properties: engagement in killing

    C = (find(L==false))';                             %check again for free spots
    M2cells= int32(C(randperm(numel(C),M2cellNo)));    %randomly allocate positions for M1s
    M2prop.Kcap = uint8(repmat(M2kmax,1,M2cellNo));    %add properties: max killing capacity
    M2prop.engaged = uint8(zeros(1,M2cellNo));         %add properties: engagement in killing
    
    if cnst.video
        figure()
    end
        
else
    %Use existing system to grow furhter
    cellfun(@(f) evalin('caller',[f ' = mySystem.TU.' f ';']), fieldnames(mySystem.TU));
    cellfun(@(f) evalin('caller',[f ' = mySystem.F.' f ';']), fieldnames(mySystem.F));
    cellfun(@(f) evalin('caller',[f ' = mySystem.M1.' f ';']), fieldnames(mySystem.M1)); 
    cellfun(@(f) evalin('caller',[f ' = mySystem.M2.' f ';']), fieldnames(mySystem.M2));
    cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));
end

   
L(TUcells)=true;                            %place tumor cells on grid
Lt = false(size(L));                        %reset tumor grid
Lt(TUcells) = true;                         %update tumor grid

L([TUcells,Fcells])=true;                   %add cells to grid
Lf = false(size(L));                        %reset fibroblast grid
Lf(Fcells) = true;                          %update fibroblast grid

L([TUcells, Fcells, M1cells])=true;         %add cells to grid
L([TUcells, Fcells, M1cells, M2cells])=true;%add cells to grid
    

%Set auxiliary variables:
nh.aux = int32([-N-1 -N -N+1 -1 1 N-1 N N+1])'; % indices to heighborhood
nh.Pms = perms(uint8(1:8))';                % permutations of adjacent positions
nh.nP = size(nh.Pms,2);                     % number of possible permutations
fcount = 0;                                 %frame counter for video
finalImage = [];                            %set empty resulting image

%ITERATIONS THROUGH NUMBER OF TIME STEPS
for i = 1:cnst.nSteps                     %iterate through time steps
    
    if (sum(M2cells)~=0 && M2TUadd > 0) %Check wether there are M2Cells present and if they are tumor promoting
        TUpprol = mySystem.params.TUpprol+M2TUadd;
        TUpprolres = mySystem.params.TUpprolres+M2TUadd;
    else %if not tumor promoting M2cells present
        TUpprol = mySystem.params.TUpprol;
        TUpprolres = mySystem.params.TUpprolres;
    end
    
    %Create fibroblast chemotaxis maps
    ChtaxMapF = double(bwdist(Lf,'euclidean'));
    
    %Tumor cell round
    L([TUcells,M1cells,M2cells,Fcells]) = true;      %ensure that all cells are on the grid
    [TUcells,TUprop] = shuffleCells(TUcells,TUprop); %shuffle tumor cells
    [L, TUcells, TUprop, TUdeath] = TU_go_grow_die(L, nh, TUcells, TUprop, TUpprol, TUpprolres, TUpmig, TUpmigres, TUpdeath, TUps, TUpres, TUrwalk, ChtaxMapF, TUpmaxres); %Update grid with growing tumor
    Lt = false(size(L));                    %reset tumor grid
    Lt(TUcells) = true;                     %update tumor grid
    %Create tumor chmotaxis map
    ChtaxMapTU = double(bwdist(Lt,'euclidean'));
     
    %Fibroblast round
    [Fcells,Fprop] = shuffleCells(Fcells,Fprop);%shuffle fibroblasts
    [L, Fcells, Fprop] = F_go_grow_die(L, nh, Fcells, Fprop, Fpprol, Fpmig, Fpdeath, Frwalk, ChtaxMapTU);
    Lf = false(size(L));                                            %reset fibroblast grid
    Lf(Fcells) = true;                                              %update fibroblast grid
    
    %M1 & M2 cell rounds
    
    [L,M1cells,M1prop,TUcells,TUprop,Lt] = mCellRound(L,Lt,M1cells,TUcells,M1prop,TUprop,M1influxProb,M1influxRate,M1kmax,M1speed,M1pmig,M1rwalk,M1pkill,M1pdeath,M1engagementDuration,ChtaxMapTU,nh,M2cells,Fcells);
    L([TUcells,M1cells,M2cells,Fcells]) = true;      %ensure that all cells are on the grid

    [L,M2cells,M2prop,TUcells,TUprop,Lt] = mCellRound(L,Lt,M2cells,TUcells,M2prop,TUprop,M2influxProb,M2influxRate,M2kmax,M2speed,M2pmig,M2rwalk,M2pkill,M2pdeath,M2engagementDuration,ChtaxMapTU,nh,M1cells,Fcells);
    L([TUcells,M1cells,M2cells,Fcells]) = true;      %ensure that all cells are on the grid
    
    mySystem=updateSystem(mySystem,TUcells,Fcells,M1cells,M2cells,TUprop,M1prop,M2prop,Fprop,L,Lt,i);
    %disp(['finished iteration ',num2str(i)]);
    
    TUcellNo(i+1,1)= length(TUcells(~TUprop.isRes)); %TUcellNo of non-resistant cells
    TUcellNo(i+1,2) = length(TUcells(TUprop.isRes)); %TUcellNo of resistant cells
    FcellNo(i+1) = length(Fcells);
    
%DRAW IF OPTION IS SELECTED
fcount=fcount+1;
if cnst.video
    visualizeSystem(TUprop,mySystem);
    drawnow, currFrame = getframe(gca);
    finalImage{fcount} = currFrame.cdata;
else
    finalImage{fcount} = false;    
end
end 
    

end

