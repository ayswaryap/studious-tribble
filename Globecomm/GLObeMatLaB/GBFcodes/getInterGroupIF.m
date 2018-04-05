function SimParams = getInterGroupIF(SimParams,sBeamM,cvxInnerM)

for iGroup = 1:SimParams.nGroups
    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
        
        groupIF = 0;
        for jGroup = 1:SimParams.nGroups
            if jGroup ~= iGroup
                for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                    xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                    groupIF = abs(SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex))^2 + groupIF;
                end
            end
        end
    end
    
    SimParams.interIF(cUserIndex,1) = groupIF;
end

end


