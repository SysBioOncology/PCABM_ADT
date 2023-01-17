function [L, Fcells, Fprop] = F_go_grow_die(L, nh, Fcells, Fprop, Fpprol, Fpmig, Fpdeath, Frwalk, ChtaxMap)

%CREATE MASKS FOR ADJACENT POSITIONS
m=getAdjacent(L,Fcells,nh);

%DETERMINE WHICH CELLS WILL PERFORM WHICH ACTIONS
% P, D and Mi are muFally exclusive; Ps and De are dependent on P
P = m.randI <= Fpprol;                                              % GROW: proliferation
D = (m.randI <= (Fpdeath+Fpprol)) & (m.randI > Fpprol);                 % DIE: spontaneous death can happen anytime
Mi = (m.randI <= (Fpdeath+Fpprol+Fpmig)) &  (m.randI > (Fpdeath+Fpprol)); % GO: migrate if no grow and no die

De = P & (Fprop.Pcap(m.indxF) == 0); % proliferation capacity exhaution -> Die
del = D | De; % cells to delete
act = find((P | Mi) & ~del); % indices to the cells that will perform action

%ADD/REMOVE ACTING CELLS TO GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % use spot with highest chemo value
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Frwalk) * chemo + Frwalk * rand(size(chemo));
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; % add new cell to grid
            if P(currID) % proliferation
                Fcells = [Fcells uint32(ngh2(indO))]; % add new cell to stack
                Fprop.Pcap(m.indxF(currID)) = Fprop.Pcap(m.indxF(currID))-1; % decrease remaining prol. cap.
                Fprop.Pcap = [Fprop.Pcap, Fprop.Pcap(m.indxF(currID))]; % update property vector for Pmax
            else % migration
                L(Fcells(m.indxF(currID))) = false; %freeing spot
                Fcells(m.indxF(currID)) = uint32(ngh2(indO));
            end
        end
    end
end

if ~isempty(del) % updating Fibroblast cell death
    L(Fcells(m.indxF(del))) = false;     % remove Fibroblast from grid
    Fcells(m.indxF(del)) = [];            % remove from stack
    Fprop.Pcap(m.indxF(del)) = [];      % remove Kmax property
end
 
end
