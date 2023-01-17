function [TUcells, TUprop, Mcells, Mprop, L, Lt] = M_promote(TUcells, TUprop, Mcells, Mprop, L, Lt,Mpprom,nh,ChtaxMap,engagementDuration)

m = getAdjacent(L,TUcells,nh);

% pre-select macrophages that may be close enough to the tumor
candidates = ChtaxMap(Mcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to promote
    K = candidates & (Mprop.engaged==0) & (Mprop.Kcap>0) & (rand(1,length(Mcells))<Mpprom);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to promote
    targetIDs = int32(zeros(1,0)); % preallocate
    promoterIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell promotion, same random order as before
    St = bsxfun(@plus,Mcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % iterate through all macrophages and look at their neighborhood
    for jj = 1:size(St,2) 
        neighbPosit = St(randperm(length(nh.aux)),jj);
        instaprom = ismember(neighbPosit(:),TUcells(:));
        possibleTargets = neighbPosit(instaprom)'; % possible targets
        targetIDs=[targetIDs,int32(possibleTargets)]; %target cells to promote
        promoterIDs = [promoterIDs, Mcells(actK(jj))]; % add promoter ID to stack
    end

    % find indices to promoted cells and promoters. 
    auxPromTU = ismember(TUcells,targetIDs); % which tumor cells are promoted
    auxPromM = ismember(Mcells,promoterIDs); % which Mmmune cells do kill

    if sum(auxPromTU)>0                 % if promoting happens, then check if promoted cells have spave & let them proliferate
        for iloop=1:numel(auxPromTU);    % for cells that can be promoted
            currID = auxPromTU(iloop);   % number within stack of currently acting cell
            ngh = m.S(:,m.indxF(currID)); % cells neighborhood
            ngh2 = ngh(ngh>0);  % erasing places that were occupied
            indO = find(~L(ngh2),1,'first');  % selecting free spot
            if ~isempty(indO)              % if there is still a free spot
                L(ngh2(indO)) = true;       % add cell to grid
                newCell = uint32(ngh2(indO)); % find place for new cell
                TUcells = [TUcells, newCell(1)]; % add new cell to stack
                TUprop.isStem = [TUprop.isStem, false];
                TUprop.Pcap = [TUprop.Pcap, TUprop.Pcap(m.indxF(currID))-1];
                if ~TUprop.isStem(m.indxF(currID)) % reduce proliferation capacity
                    TUprop.Pcap(m.indxF(currID)) = TUprop.Pcap(m.indxF(currID))-1;
                end
                if  TUprop.isRes(m.indxF(currID))
                    TUprop.isRes = [TUprop.isRes, true];
                else
                    TUprop.isRes = [TUprop.isRes, false];
                end
                
            end
        end
    end
    Mprop.engaged(auxPromM) = engagementDuration; % promoters are engaged
    %Mprop.Kcap(auxPromM) = Mprop.Kcap(auxPromM)-1; % exhaust promoters
    L(TUcells) = true;  % FIRST make sure all cels are on L grid
    Lt(TUcells) = true;  % ... and on Lt grid
    end % end actual killing filter
end


end