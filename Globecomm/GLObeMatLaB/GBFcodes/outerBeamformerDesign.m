
function SimParams = outerBeamformerDesign(SimParams)

% * * Depending upon the choice of outer BF we design the outer BF * *

gCovMatrix = cell(SimParams.nGroups,1);
covMatrix = zeros(SimParams.nTransmit);
for iGroup = 1:SimParams.nGroups    
    gCovMatrix{iGroup,1} = zeros(SimParams.nTransmit);
    for iRealization = 1:SimParams.nRealizations
        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
            cH = getRingChannel(SimParams,iGroup,iUser,'RING');
            covMatrix = cH' * cH + covMatrix;
            gCovMatrix{iGroup,1} = gCovMatrix{iGroup,1} + cH' * cH;
        end
    end    
end

switch upper(SimParams.statBeamType)
    
    case {'EIG'} % Statistical information
        
        [U,~,~] = svd(covMatrix);
        for iGroup = 1:SimParams.nGroups % When users have full access R
            beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
            SimParams.groupInfo(iGroup).statBeams = U(:,beamIndices);
        end
        
    case 'EIG_G'
        for iGroup = 1:SimParams.nGroups % When group specific processing R_g
            [U,~,~] = svd(gCovMatrix{iGroup,1});
            SimParams.groupInfo(iGroup).statBeams = U(:,1:SimParams.gStatBeams);
        end
                
    case 'EYE' % Digital beamformer
        
        for iGroup = 1:SimParams.nGroups %No outer beamformer
            SimParams.limitToGroupBeamsOnly = 0;
            SimParams.groupInfo(iGroup).activeAntennas = ones(size(SimParams.nTransmit));
            SimParams.groupInfo(iGroup).statBeams = eye(SimParams.nTransmit);
        end
        SimParams.gStatBeams = SimParams.nTransmit;
        
    case 'UDFT' % Based on user location information (Group specific covariance matrix)
        
        for iGroup = 1:SimParams.nGroups
            covMatrix = zeros(SimParams.nTransmit);
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cH = getRingChannel(SimParams,iGroup,iUser,'AWGN');
                covMatrix = cH' * cH + covMatrix;
            end
            SimParams.groupInfo(iGroup).covMatrix = covMatrix;
        end
        
        for iGroup = 1:SimParams.nGroups
            [U,~,~] = svd(SimParams.groupInfo(iGroup).covMatrix);
            SimParams.groupInfo(iGroup).statBeams = U(:,1:SimParams.gStatBeams);
        end
        
    case 'BDFT' % Based on max and min of user location [uniformly sampled]
        
        oversampFactor = 2^0;
        U = zeros(SimParams.nTransmit,SimParams.gStatBeams * oversampFactor);
        couplingMatrix = kron(eye(SimParams.gStatBeams),ones(oversampFactor,1));
        for iGroup = 1:SimParams.nGroups
            thetaBegin = min(SimParams.groupUserBaseAngles{iGroup,1});
            thetaEnd = max(SimParams.groupUserBaseAngles{iGroup,1});
            beamDirections = linspace(thetaBegin,thetaEnd,oversampFactor * SimParams.gStatBeams + 1);
            beamDirections = beamDirections(1:end-1);
            for iBeam = 1:length(beamDirections)
                U(:,iBeam) = exp(-sqrt(-1) * pi * sind(beamDirections(iBeam)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit)';
            end
            SimParams.groupInfo(iGroup).statBeams = U * couplingMatrix;
            SimParams.groupInfo(iGroup).beamDirections = beamDirections;
        end
        
        covMatrix = zeros(SimParams.nTransmit);
        for iGroup = 1:SimParams.nGroups
            covMatrix = SimParams.groupInfo(iGroup).statBeams * SimParams.groupInfo(iGroup).statBeams' + covMatrix;            
        end
        
        [U,~,~] = svd(covMatrix);
        for iGroup = 1:SimParams.nGroups
            beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
            SimParams.groupInfo(iGroup).statBeams = U(:,beamIndices);
        end
        
    case 'DFT' % Based on DFTMTX (subsampling)
        
        dftMatrix = dftmtx(SimParams.nTransmit);        
        subIndices = round(linspace(1,SimParams.nTransmit,SimParams.totStatBeams));
        U = dftMatrix(:,subIndices);
        for iGroup = 1:SimParams.nGroups
            beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
            SimParams.groupInfo(iGroup).statBeams = U(:,beamIndices);
        end
                        
    case 'GDFT' % Divides entire space uniformly over [-90, 90] depending on total number of beams
        
        oversampFactor = 2^0;
        totBeamDirections = reshape(linspace(-90,90,SimParams.totStatBeams),[],SimParams.nGroups).';
        U = zeros(SimParams.nTransmit,SimParams.gStatBeams * oversampFactor);
        couplingMatrix = kron(eye(SimParams.gStatBeams),ones(oversampFactor,1));
        for iGroup = 1:SimParams.nGroups
            beamDirections = totBeamDirections(iGroup,:);
            for iBeam = 1:length(beamDirections)
                U(:,iBeam) = exp(-sqrt(-1) * pi * sind(beamDirections(iBeam)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit)';
            end
            SimParams.groupInfo(iGroup).statBeams = U * couplingMatrix;
            SimParams.groupInfo(iGroup).beamDirections = beamDirections;
        end
        
        covMatrix = zeros(SimParams.nTransmit);
        for iGroup = 1:SimParams.nGroups
            covMatrix = SimParams.groupInfo(iGroup).statBeams * SimParams.groupInfo(iGroup).statBeams' + covMatrix;            
        end
        
        [U,~,~] = svd(covMatrix);
        for iGroup = 1:SimParams.nGroups
            beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
            SimParams.groupInfo(iGroup).statBeams = U(:,beamIndices);
        end
        
    case 'GREEDY'
        
        dftMatrix = dftmtx(SimParams.nTransmit);
        beamProjGains = zeros(SimParams.nTransmit,1);
        for iBeam = 1:SimParams.nTransmit
            beamProjGains(iBeam,1) = real(trace(dftMatrix(:,iBeam)' * covMatrix * dftMatrix(:,iBeam)));
        end
        [~,iInd] = sort(beamProjGains,'descend');

        for iGroup = 1:SimParams.nGroups
            beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
            SimParams.groupInfo(iGroup).statBeams = dftMatrix(:,iInd(beamIndices));
        end
        
    case 'GREEDY_G'
        dftMatrix = dftmtx(SimParams.nTransmit);
        beamProjGains = zeros(SimParams.nTransmit,1);
        for iGroup = 1:SimParams.nGroups
            for iBeam = 1:SimParams.nTransmit
                beamProjGains(iBeam,1) = real(trace(dftMatrix(:,iBeam)' * gCovMatrix{iGroup,1} * dftMatrix(:,iBeam)));
            end
            [~,iInd] = sort(beamProjGains,'descend');
            SimParams.groupInfo(iGroup).statBeams = dftMatrix(:,iInd(1:SimParams.gStatBeams));
        end
        
        
%         for iGroup = 1:SimParams.nGroups
%             beamIndices = ((iGroup - 1) * SimParams.gStatBeams + 1) : (iGroup * SimParams.gStatBeams);
%             SimParams.groupInfo(iGroup).statBeams = dftMatrix(:,iInd(beamIndices));
%         end

        
    case 'OTHER'
        if strcmp(opcmode,'greedy1')
            [U eigval] = eig(R);
            k = 1;
            selected = 1;
            while selected <= S
                [~,n] = size(U);
                clear bestforuserk;
                for u=1:n
                    bestforuserk(u) = U(:,u)'*R_k(:,:,k)*U(:,u);
                end
                [~,ind] = max(bestforuserk);
                B(:,selected) = U(:,ind);
                U(:,ind) = [];
                if k < K
                    k = k+1;
                else
                    k = 1;
                end
                selected = selected + 1;
            end
        end
        
        if strcmp(opcmode, 'greedy2')
            [U eigval] = eig(R);
            for u=1:M
                for k=1:K
                    bestcorr(k,u) = U(:,u)'*R_k(:,:,k)*U(:,u);
                end
            end
            selected = 1;
            while selected <= S
                [maxperrow,row] = max(bestcorr);
                [maxcorr,column] = max(maxperrow);
                B(:,selected) = U(:,column);
                bestcorr(:,column) = [];
                U(:,column) = [];
                selected = selected + 1;
            end
        end
        
end

%Depending on the choice of outer BF we design the B accordingly

if ~strcmpi(SimParams.statBeamType,'EYE')
    SimParams.sBeamM = [];
    activeBeamsPerGroup = zeros(1,SimParams.totStatBeams);
    activeBeamsPerGroup(1:SimParams.gStatBeams) = 1;
    for iGroup = 1:SimParams.nGroups
        SimParams.groupInfo(iGroup).activeBeams = circshift(activeBeamsPerGroup, SimParams.gStatBeams * (iGroup - 1), 2);
        SimParams.sBeamM = [SimParams.sBeamM, SimParams.groupInfo(iGroup).statBeams];
    end
else
    fprintf('Outer precoder is not defined !! \n');
    SimParams.sBeamM = eye(SimParams.nTransmit);
end

end