function [TUcells, TUprop, Mcells, Mprop, L, Lt] = M_kill(TUcells, TUprop, Mcells, Mprop, L, Lt,Mpkill,nh,ChtaxMap,engagementDuration)

% pre-select macrophages that may be close enough to the tumor
candidates = ChtaxMap(Mcells)<=1;
if sum(candidates(:)) % if there are candidates
    % select cells that are going to kill
    K = candidates & (Mprop.engaged==0) & (Mprop.Kcap>0) & (rand(1,length(Mcells))<Mpkill);
    actK = find(K); % cell indices
    if ~isempty(actK) % if there is a cell that is going to kill
    targetIDs = int32(zeros(1,0)); % preallocate
    killerIDs = int32(zeros(1,0)); % preallocate
    % start tumor cell killing, same random order as before
    St = bsxfun(@plus,Mcells(actK),nh.aux(nh.Pms(:,randi(nh.nP,1,length(actK)))));
    % iterate through all macrophages and look at their neighborhood
    for jj = 1:size(St,2) 
        neighbPosit = St(randperm(length(nh.aux)),jj);
        instakill = ismember(neighbPosit(:),TUcells(:));
        % if the cell encounters another cell to kill
        if sum(instakill) > 0 
            % if more than 1 possible targets then use the first one
            possibleTargets = neighbPosit(instakill); % possible targets
            killwhat = int32(possibleTargets(1)); % kill only the first candidate   
            targetIDs = [targetIDs, killwhat]; % add target ID to stack
            killerIDs = [killerIDs, Mcells(actK(jj))]; % add killer ID to stack
        end
    end

    % find indices to killed cell and killer cell. If the unlikely case
    % happens that one tumor cell is simultaneously killed by two macrophages,
    % then both will be exhausted
    auxKillTU = ismember(TUcells,targetIDs); % which tumor cells are killed
    auxKillM = ismember(Mcells,killerIDs); % which Mmmune cells do kill

    if sum(auxKillTU)>0                 % if killing happens, then update  
        L(TUcells(auxKillTU)) = false;  % FIRST remove from L grid
        Lt(TUcells(auxKillTU)) = false;  % ... and remove from Lt grid
        TUcells(auxKillTU) = [];            % remove from stack
        TUprop.isStem(auxKillTU) = [];      % remove stemness property
        TUprop.isRes(auxKillTU) =[];        % remove resistant property
        TUprop.Pcap(auxKillTU) = [];        % remove Pmax property
        Mprop.Kcap(auxKillM) = Mprop.Kcap(auxKillM)-1; % exhaust killers
        Mprop.engaged(auxKillM) = engagementDuration; % killers are engaged
    end

    end % end actual killing filter
end % end candidate filter


end