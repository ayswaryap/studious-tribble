function SimParams = getIntraGroupIF(SimParams,sBeamM,cvxInnerM)

for iGroup = 1:SimParams.nGroups
    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
        
        groupIF = sqrt(SimParams.N0);
        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
            if jUser ~= iUser
                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                groupIF = groupIF + abs(SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex))^2;
            end
        end
        
        SimParams.intraIF(cUserIndex,1) = groupIF;
    end
end

end
