function [L, Mcells, Mprop] =  M_go_die(Mcells, Mprop, Mpmig, Mpdeath, Mrwalk, ChtaxMap, L, nh)
    
%CREATE MASKS FOR ADJACENT POSITIONS
m = getAdjacent(L,Mcells,nh);

%DETERMINE WHICH CELLS PERFORM ACTIONS
% D and Mi are mutually exclusive
D = m.randI <= (Mpdeath);                                        % DIE: spontaneous death can happen anytime
Mi = (m.randI <= (Mpdeath+Mpmig)) &  (m.randI > (Mpdeath));      % GO: migrate if no grow and no die

del = D; % cells to delete
act = find(Mi & ~del); % indices to the cells that will perform action

%MOVE CELLS ON GRID
for iloop = 1:numel(act) % only for those that will do anything
    currID = act(iloop); % number within stack of currently acting cell
    ngh = m.S(:,m.indxF(currID)); % cells neighborhood
    ngh2 = ngh(ngh>0); % erasing places that were previously occupied
    indOL = find(~L(ngh2)); %selecting all free spots  
    chemo = ChtaxMap(ngh2(:)); % extract chemotaxis value at neighbors
    chemo = chemo(~L(ngh2)); % block occupied spots   
    if ~isempty(chemo) % use spot with highest chemo value
        chemo = chemo/max(chemo(:)); % normalize
        chemo = (1-Mrwalk) * chemo + Mrwalk * rand(size(chemo));
        [~,cid] = min(chemo); % lowest distance 
        indO = indOL(cid(1));   
        if ~isempty(indO) %if there is still a free spot
            L(ngh2(indO)) = true; % add new cell to grid
            L(Mcells(m.indxF(currID))) = false; %freeing spot
            Mcells(m.indxF(currID)) = uint32(ngh2(indO));
        end
    end
end

if ~isempty(del) % updating macrophage cell death
    L(Mcells(m.indxF(del))) = false;    % remove macrophage cell from grid
    Mcells(m.indxF(del)) = [];          % remove from stack
    Mprop.Kcap(m.indxF(del)) = [];      % remove Kmax property
    Mprop.engaged(m.indxF(del)) = [];   % remove engagement property
end
 
end
