function [L, TUcells, TUprop, TUdeath] = TU_go_grow_die(L, nh, TUcells, TUprop, TUpprol, TUpprolres, TUpmig, TUpmigres, TUpdeath, TUps, TUpres, TUrwalk, ChtaxMap, TUpmaxres)

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,TUcells,nh);

%DETERMINE WHICH CELLS WILL PERFORM WHICH ACTIONS
% P, D and Mi are mutually exclusive; Ps and De are dependent on P

if sum(TUprop.isRes) > 0                                   %Check if there are any resistant cells
    %Create seperate vectors for normal and resistant cells
    %Vectors are the same, but there is a 1 in case of non/resistant cell (1
    %means that no action will be performed as 1 is never smaller than the
    %set probability.
    isRes = TUprop.isRes(m.indxF);                        %Resistant cells with free space
    randiNorm = m.randI; randiNorm(isRes) = 1;
    randiRes = m.randI; randiRes(~isRes) = 1;
      
    % GROW: proliferation of Normal and Resistant cells
    Pnorm = randiNorm <= TUpprol;                                            
    Pres = randiRes <= TUpprolres;
    Pres = Pnorm | Pres;                                           
    
    % DIE: spontaneous death can happen anytime
    Dnorm = (randiNorm <= (TUpdeath+TUpprol)) & (randiNorm > TUpprol);  
    Dres = (randiRes <= (TUpdeath+TUpprolres)) & (randiRes > TUpprolres); 
    D = Dnorm | Dres;
    
    % GO: migrate if no grow and no die
    Minorm = (randiNorm <= (TUpdeath+TUpprol+TUpmig)) &  (m.randI > (TUpdeath+TUpprol));
    Mires = (randiRes <= (TUpdeath+TUpprolres+TUpmigres)) &  (m.randI > (TUpdeath+TUpprolres));
    Mi = Minorm | Mires;
    
    Ps = Pres & rand(1,m.nC) <= TUps & TUprop.isStem(m.indxF);               % symmetric division
    Pm = Pres & ~Ps & rand(1,m.nC) <= TUpres & ~TUprop.isRes(m.indxF);                            % cell generates resistant cell upon division (does not happen for symmetric division)
else
    Pres = m.randI <= TUpprol;                                                    % GROW: proliferation
    D = (m.randI <= (TUpdeath+TUpprol)) & (m.randI > TUpprol);                 % DIE: spontaneous death can happen anytime
    Mi = (m.randI <= (TUpdeath+TUpprol+TUpmig)) &  (m.randI > (TUpdeath+TUpprol)); % GO: migrate if no grow and no die
    
    Ps = Pres & rand(1,m.nC) <= TUps & TUprop.isStem(m.indxF);         % symmetric division
    Pm = Pres & rand(1,m.nC) <= TUpres & ~Ps;                          % cell creates resistant cell upon devision (does not happen for symmetric division)
end

De = Pres & (TUprop.Pcap(m.indxF) == 0); %& (TUprop.isRes(m.indxF) == 0);                       % proliferation capacity exhaustion -> Die
del = D | De;                                               % find dead / dying cells
act = find((Pres | Mi) & ~del);                                % live cells that will proliferate or migrate

%ADD/REMOVE ACTING CELLS TO GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % if there is still a free spot
        chemo = (1-TUrwalk) * chemo + TUrwalk * rand(size(chemo));
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; % add cell to grid
            if Pres(currID) % proliferation happens
                newCell = uint32(ngh2(indO)); % find place for new cell
                TUcells = [TUcells, newCell(1)]; % add new cell to stack
                if Ps(currID) % symmetric division
                    TUprop.isStem = [TUprop.isStem, true];
                    TUprop.isRes = [TUprop.isRes, false];
                    TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))];  
                elseif Pm(currID) %resistant cell division
                    TUprop.isStem = [TUprop.isStem, false];
                    TUprop.isRes = [TUprop.isRes, true];
                    %TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))-1]; %TUprop.Pcap(m.indxF(currID))-1
                    %TUprop.Pcap = [TUprop.Pcap, Inf];
                    %TUprop.Pcap = [TUprop.Pcap, 50];
                    TUprop.Pcap = [TUprop.Pcap, TUpmaxres];
                    if ~TUprop.isStem(m.indxF(currID)) % reduce proliferation capacity
                        TUprop.Pcap(m.indxF(currID)) = TUprop.Pcap(m.indxF(currID))-1;
                    end
                else % asymmetric division
                    TUprop.isStem = [TUprop.isStem, false];
                    TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))-1]; %TUprop.Pcap(m.indxF(currID))-1
                    if  TUprop.isRes(m.indxF(currID))
                        TUprop.isRes = [TUprop.isRes, true];
                    else
                        TUprop.isRes = [TUprop.isRes, false];
                    end
                    if ~TUprop.isStem(m.indxF(currID)) % reduce proliferation capacity
                        TUprop.Pcap(m.indxF(currID)) = TUprop.Pcap(m.indxF(currID))-1;
                    end
                end
            else % migration
                L(TUcells(m.indxF(currID))) = false; % freeing spot
                TUcells(m.indxF(currID)) = uint32(ngh2(indO)); % update cell position
            end
        end
    end
end

if ~isempty(del) % remove dead tumor cells
    L(TUcells(m.indxF(del))) = false;      % remove from grid
    TUcells(m.indxF(del)) = [];            % remove from stack
    TUprop.isStem(m.indxF(del)) = [];      % remove stemness property
    TUprop.isRes(m.indxF(del)) = [];       % remove resistant property
    TUprop.Pcap(m.indxF(del)) = [];        % remove Pmax property
    
end 

TUdeath = sum(D);                        % number of dying tumor cells

end
