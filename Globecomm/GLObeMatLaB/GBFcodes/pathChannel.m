function [varargout] = pathChannel(SimParams, iGroup, iUser, varargin)

theta = 2*pi*SimParams.groupInfo(iGroup).userLocs(iUser);
Pathloss_range = [0, 20];
pathloss_range_lin = 10.^(Pathloss_range/10);
SimParams.groupInfo(iGroup).nUsers = SimParams.groupInfo(iGroup).nUsers;
SimParams.groupInfo(iGroup).angSpread = SimParams.angSpread * ones(1,SimParams.groupInfo(iGroup).nUsers);

if size(varargin,2)
    
    switch varargin{1}
        
        case 'Initialize'
            
            
            pathloss = (pathloss_range_lin(2) - pathloss_range_lin(1))*rand(SimParams.groupInfo(iGroup).nUsers,1) + pathloss_range_lin(1); % user path losses
            
            for paths = 1:SimParams.numScatterers
                
                relative_path_angle = (SimParams.groupInfo(iGroup).angSpread(iUser)/180)*pi*rand(1) - (SimParams.groupInfo(iGroup).angSpread(iUser)/360)*pi;
                thetapath = theta + relative_path_angle; % path AoD
                
                if thetapath < 0 % all angles between 0 - 2pi
                    thetapath = 2*pi + thetapath;
                elseif thetapath > 2*pi
                    thetapath = thetapath - 2*pi;
                end
                
                for nTx = 1:SimParams.nTransmit % Antenna signature vector for path p
                    antennaSignature(nTx) = exp(-1i*pi*(nTx-1)*cos(thetapath));
                end
                
                pathchan(SimParams.numScatterers,:) = antennaSignature*exp(1i*2*pi*rand(1)); % add random phase noise
                
            end
            
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1} = sqrt(pathloss(iUser)/SimParams.numScatterers)*sum(pathchan);
            
            varargout{1} = SimParams;
            
        case 'Reset'
            
            pathloss = (pathloss_range_lin(2) - pathloss_range_lin(1))*rand(SimParams.groupInfo(iGroup).nUsers,1) + pathloss_range_lin(1); % user path losses
            
            for paths = 1:SimParams.numScatterers
                
                relative_path_angle = (SimParams.groupInfo(iGroup).angSpread(iUser)/180)*pi*rand(1) - (SimParams.groupInfo(iGroup).angSpread(iUser)/360)*pi;
                thetapath = theta + relative_path_angle; % path AoD
                
                if thetapath < 0 % all angles between 0 - 2pi
                    thetapath = 2*pi + thetapath;
                elseif thetapath > 2*pi
                    thetapath = thetapath - 2*pi;
                end
                
                for nTx = 1:SimParams.nTransmit % Antenna signature vector for path p
                    antennaSignature(nTx) = exp(-1i*pi*(nTx-1)*cos(thetapath));
                end
                
                pathchan(SimParams.numScatterers,:) = antennaSignature*exp(1i*2*pi*rand(1)); % add random phase noise
                
            end
            
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1} = sqrt(pathloss(iUser)/SimParams.numScatterers)*sum(pathchan);
            
            
            for i = 1:length(SimParams.nGroups)
                varargout{i} = cell2mat(SimParams.groupInfo(iGroup).mimoChannel(i));
            end
            
            
            
    end
    
    
else
    
    display('There is a problem ');
    
end

end

