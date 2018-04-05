function ZF = initializeSCAPoints(SimParams,cvxType)

sBeamM = SimParams.sBeamM;

switch cvxType
    
    case 'CVX'
        
        H = zeros(SimParams.nUsers,SimParams.nTransmit);
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser);
            end
        end
        
        effChannel = H * sBeamM;
        ZF = effChannel' / (effChannel * effChannel');
        
        zfPower = diag(1./sqrt(diag(ZF' * ZF)));
        %         wfPower = waterfill(SimParams.txPower,reshape(diag(ZF' * ZF),1,[]));
        wfPower = SimParams.txPower;
        ZF = ZF * sqrt(diag(wfPower)) * zfPower;
        eZF = sBeamM * ZF;
        ZF = ZF * sqrt(SimParams.txPower) / sqrt(trace(eZF * eZF'));
        
    case 'CVXG'
        
        ZF = zeros(SimParams.gStatBeams,SimParams.nUsers);
        effChannel = zeros(SimParams.nUsers,SimParams.gStatBeams);
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                effChannel(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                ZF(:,cUser) = effChannel(cUser,:)' / norm(effChannel(cUser,:));
            end
        end
        
        ZF = ZF * (sqrt(SimParams.txPower / sqrt(real(trace(ZF * ZF')))));
        
    case 'KKT'
        
        ZF = zeros(SimParams.gStatBeams,SimParams.nUsers);
        effChannel = zeros(SimParams.nUsers,SimParams.gStatBeams);
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                effChannel(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                ZF(:,cUser) = effChannel(cUser,:)' / norm(effChannel(cUser,:));
            end
        end
        
        ZF = ZF * (sqrt(SimParams.txPower / sqrt(real(trace(ZF * ZF')))));
        
        
    case 'CVXSOC'
        
        H = zeros(SimParams.nUsers,SimParams.nTransmit);
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser);
            end
        end
        
        effChannel = H * sBeamM; % * zeroMatrix;
        %         ZF = effChannel' / (effChannel * effChannel');
        
        ZF = effChannel' * diag(1 ./ sqrt(diag(effChannel * effChannel')));
        
        zfPower = diag(1./sqrt(diag(ZF' * ZF)));
        wfPower = waterfill(SimParams.txPower,reshape(diag(ZF' * ZF),1,[]));
        ZF = ZF * sqrt(diag(wfPower)) * zfPower;
        
end

end
