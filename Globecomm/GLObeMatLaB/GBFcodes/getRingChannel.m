function channelH = getRingChannel(SimParams, iGroup, iUser, chnType, varargin)

cUser = SimParams.groupUserIndices{iGroup,1}(1,iUser);
userLoc = SimParams.groupInfo(iGroup).userLocs(1,iUser);

switch upper(chnType)
    case 'AWGN'
        nTaps = 1;
        channelH = zeros(nTaps,SimParams.nTransmit);
        for iTap = 1:nTaps
            channelH(iTap,:) = exp(-sqrt(-1) * pi * sind(userLoc) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
        end
    case 'LOS'
        nTaps = 1;
        cplxGains = ones(1,nTaps);
        channelH = zeros(nTaps,SimParams.nTransmit);
        for iTap = 1:nTaps
            channelH(iTap,:) = exp(-sqrt(-1) * pi * sind(userLoc) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
        end
        channelH = cplxGains * channelH * sqrt(SimParams.pathLoss(1,cUser) / nTaps);
    case 'RING'
        nTaps = SimParams.numScatterers;
        cplxGains = exp(sqrt(-1) * rand(1,nTaps) * 2 * pi);
        userLoc = userLoc + (2 * rand(1,nTaps) - 1) * SimParams.angSpread;
        channelH = zeros(nTaps,SimParams.nTransmit);
        for iTap = 1:nTaps
            channelH(iTap,:) = exp(-sqrt(-1) * pi * sind(userLoc(1,iTap)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
        end
        channelH = cplxGains * channelH * sqrt(SimParams.pathLoss(1,cUser) / nTaps);
        
    case 'UCA'
        nTaps = SimParams.numScatterers;
        cplxGains = exp(sqrt(-1) * rand(1,nTaps) * 2 * pi);
        userElevAngle = SimParams.groupInfo(iGroup).elevationAngle(1,iUser);
        userLoc = userLoc + (2 * rand(1,nTaps) - 1) * SimParams.angSpread;
        channelH = zeros(nTaps,SimParams.nTransmit);
        for iTap = 1:nTaps
            channelH(iTap,:) = exp(-sqrt(-1) * pi * sind(userElevAngle) * cosd(userLoc(1,iTap)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
        end
        channelH = cplxGains * channelH * sqrt(SimParams.pathLoss(1,cUser) / nTaps);
                
    case 'TEMP'
        
        if strcmp(varargin,'Init')
            nTaps = SimParams.numScatterers;
            cplxGains = complex(randn(1,nTaps),randn(1,nTaps));
            userLoc = userLoc + (2 * rand(1,nTaps) - 1) * SimParams.angSpread;
            cplxGains = cplxGains / norm(cplxGains);
            channelHi = zeros(nTaps,SimParams.nTransmit);
            for iTap = 1:nTaps
                channelHi(iTap,:) = exp(-sqrt(-1) * pi * sind(userLoc(1,iTap)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
            end
            channelH = channelHi;
            
        elseif strcmp(varargin,'Reset')
            nTaps = SimParams.numScatterers;
            cplxGains = complex(randn(1,nTaps),randn(1,nTaps));
            userLoc = userLoc + (2 * rand(1,nTaps) - 1) * SimParams.angSpread;
            cplxGains = cplxGains / norm(cplxGains);
            channelHi = zeros(nTaps,SimParams.nTransmit);
            for iTap = 1:nTaps
                channelHi(iTap,:) = exp(-sqrt(-1) * pi * sind(userLoc(1,iTap)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
            end
            channelH = sum(channelHi);
            
        end
        
end

channelH = channelH * sqrt(SimParams.nTransmit);

end