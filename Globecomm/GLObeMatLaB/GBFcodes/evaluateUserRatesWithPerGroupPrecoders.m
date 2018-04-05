function SimParams = evaluateUserRatesWithPerGroupPrecoders(SimParams,cvxInnerM)

SimParams.totUserRateE = zeros(SimParams.nUsers,1);
SimParams.totUserBeta = zeros(SimParams.nUsers,1);


for iGroup = 1:SimParams.nGroups
    
    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
        
        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
        totalIF = sqrt(SimParams.N0);
        
        %Interference from my own group: Intragroup
        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
            if jUser ~= iUser
                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,xUserIndex)];
            end
        end
        
        %Interference from my neighbouring group: Intergroup
        for jGroup = 1:SimParams.nGroups
            if jGroup ~= iGroup
                for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                    xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                    totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(jGroup).statBeams * cvxInnerM(:,xUserIndex)];
                end
            end
        end
        
        userBeta = totalIF * totalIF';
        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
        userGamma = abs(effChannel * cvxInnerM(:,cUserIndex))^2 / userBeta;
        SimParams.totUserRateE(cUserIndex,1) = log2(1 + userGamma);
        SimParams.totUserBeta(cUserIndex,1) = userBeta;
        
    end
    
end

end


