function [varargout] = getChannel(SimParams, iGroup, iUser, varargin)

if size(varargin,2)
    
    switch varargin{1}
        case 'Initialize'            
            
            uTheta = SimParams.groupInfo(iGroup).userLocs(iUser) * pi / 180;
            [uLoc.x, uLoc.y, uLoc.z] = sph2cart(uTheta,0,SimParams.userRingRadius);
            
            sLoc = struct('x',0,'y',0,'z',0);
            scatPositions = SimParams.angSpread * (rand(1,SimParams.numScatterers) - 0.5) * pi / 180 + uTheta;
            for iScatterer = 1:SimParams.numScatterers
                [sLoc(iScatterer).x, sLoc(iScatterer).y, sLoc(iScatterer).z] = sph2cart(scatPositions(1,iScatterer),0,(SimParams.userRingRadius + (rand - 0.5) * 5));
            end
            sGains = complex(randn(1,SimParams.numScatterers),randn(1,SimParams.numScatterers));
            
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1} = ...
                phased.ScatteringMIMOChannel(...
                'TransmitArray',SimParams.UCA,...
                'ReceiveArray',SimParams.rxArray,...
                'PropagationSpeed',physconst('LightSpeed'),...
                'CarrierFrequency',SimParams.cFrequency,...
                'SpecifyAtmosphere',false,...
                'SampleRate',SimParams.sFrequency,...
                'TransmitArrayPosition',[0;0;0],...
                'ReceiveArrayPosition',reshape(struct2array(uLoc),3,[]),...
                'NumScatterers',SimParams.numScatterers,...
                'ScattererSpecificationSource','Property',...
                'ScattererPosition',reshape(struct2array(sLoc),3,[]),...
                'ScattererCoefficient',sGains,...
                'SimulateDirectPath',true,'ChannelResponseOutputPort',true,...
                'ReceiveArrayMotionSource','Input port');
            
            SimParams.groupInfo(iGroup).userPlatform{iUser,1} = phased.Platform('MotionModel','Velocity','InitialPosition',...
                SimParams.groupInfo(iGroup).mimoChannel{iUser,1}.ReceiveArrayPosition,'Velocity', 2 * rand(3,1) - 1);
        
        case 'Reset'
            
            sLoc = struct('x',0,'y',0,'z',0);
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1}.release();
            uTheta = SimParams.groupInfo(iGroup).userLocs(iUser) * pi / 180;
            
            scatPositions = SimParams.angSpread * (rand(1,SimParams.numScatterers) - 0.5) * pi / 180 + uTheta;
            for iScatterer = 1:SimParams.numScatterers
                [sLoc(iScatterer).x, sLoc(iScatterer).y, sLoc(iScatterer).z] = sph2cart(scatPositions(1,iScatterer),0,(SimParams.userRingRadius + (rand - 0.5) * 5));
            end
            sGains = complex(randn(1,SimParams.numScatterers),randn(1,SimParams.numScatterers));
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1}.ScattererCoefficient = sGains; 
            SimParams.groupInfo(iGroup).mimoChannel{iUser,1}.ScattererPosition = reshape(struct2array(sLoc),3,[]);
    end
    
    varargout{1} = SimParams;

else    
    mimoChannel = SimParams.groupInfo(iGroup).mimoChannel{iUser};
    [rxPos,rxVel] = SimParams.groupInfo(iGroup).userPlatform{iUser,1}(1e-2);
    [~, chanMatrix, ~] = step(mimoChannel,ones(1,SimParams.nTransmit),rxPos,rxVel,rotz(randi(360,1,1)));
    varargout{1} = sum(chanMatrix,3) * db2pow(20);

    if 0
        [rxPos,rxVel] = SimParams.groupInfo(iGroup).userPlatform{iUser,1}(1e-1);
        [~, chanMatrix, ~] = step(mimoChannel,ones(1,SimParams.nTransmit),rxPos,rxVel,rotz(randi(360,1,1)));
        plot(real(sum(chanMatrix,3)));hold all;
    end
    
%     angRange = [-179 : 180];
%     SV = step(SimParams.stVector,SimParams.cFrequency,angRange);
%     plot(angRange,abs(SV' * varargout{1}));
    
%     angRange = [-179 : 180];
%     SV = step(SimParams.stVector,SimParams.cFrequency,angRange)';
%     sProj = step(mimoChannel,SV);
%     plot(angRange,abs(sProj));
end

end











% function [pathChannel] = getChannel(SimParams,baseTheta)
%
% cTheta = (2 * rand(1,SimParams.nPaths) - 1) * SimParams.angSpread; %angles
% cTheta = wrapTo180(cTheta + baseTheta);
%
% fs = 10e6;
% Ns = 0;
%
%
%
%
% Channel = phased.ScatteringMIMOChannel(...
%     'TransmitArray',SimParams.txAntennas,...
%     'ReceiveArray',SimParams.rxAntennas,...
%     'PropagationSpeed',SimParams.c,...
%     'CarrierFrequency',SimParams.frequency,...
%     'SpecifyAtmosphere',true,...
%     'SampleRate',fs,...
%     'TransmitArrayPosition',[0;0;0],...
%     'ReceiveArrayPosition',[cTheta;0;0],...
%     'NumScatterers',Ns,'ScattererPosition',[50 ; -50 ; 0],...
%     'SimulateDirectPath',true,'ChannelResponseOutputPort',false,...
%     'ScattererSpecificationSource','Property',...
%     'SeedSource','Property','Seed',0);
%
%
%
%
%
%  SV = step(SimParams.steerVector,SimParams.frequency,SimParams.groupInfo.baseTheta).';
%
%  pathChannel = step(Channel,SV);
%
%
%  %[x,y] = rangeangle(Channel.ScattererPosition,Channel.TransmitArrayPosition);
% %
%  plot(SimParams.groupInfo.baseTheta,db(abs(pathChannel)));
%
% % elementLocs = linspace(0,(SimParams.nTransmit - 1),SimParams.nTransmit) * 0.5;
% % pathChannel = cbfweights(elementLocs,cTheta) * diag(exp(sqrt(-1) * 2 * pi * rand(1,SimParams.nPaths)) .* 10.^(-(1:SimParams.nPaths)/20));
% % pathChannel = sum(pathChannel,2) / sqrt(sum(10.^(-(1:SimParams.nPaths)/20)));
% % pathChannel = shiftdim(pathChannel,1);
% %
% %
% % resp = step(SimParams.response,SimParams.frequency,SimParams.groupInfo.baseTheta,sv);
% %
%
% end
%
