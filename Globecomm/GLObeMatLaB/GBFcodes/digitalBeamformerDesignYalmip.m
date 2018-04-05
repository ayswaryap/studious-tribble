
function SimParams = digitalBeamformerDesignYalmip(SimParams)

if ispc
    gprintf = @fprintf;
else
    gprintf = @nprintf;
end

alphaM = 0.5;
sBeamM = SimParams.sBeamM;
precType = strsplit(SimParams.innerPrecoder,'_');

switch precType{1}
    
    case 'ZF'
        
        for iRun = 1:SimParams.nMontRuns
            
            H = zeros(SimParams.nUsers,SimParams.nTransmit);
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                    H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser);
                end
            end
            
            effChannel = H * sBeamM; % * zeroMatrix;
            ZF = effChannel' / (effChannel * effChannel');
            
            zfPower = diag(1./sqrt(diag(ZF' * ZF)));
            wfPower = SimParams.txPower;
            %             wfPower = waterfill(SimParams.txPower,reshape(diag(ZF' * ZF),1,[]));
            ZF = ZF * sqrt(diag(wfPower)) * zfPower;
            eZF = sBeamM * ZF;
            cvxInnerMP = ZF * sqrt(SimParams.txPower) / sqrt(trace(eZF * eZF'));
            SimParams.sBeamM = sBeamM;
            SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
            
            SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
            
            SimParams.groupSumRate.srate(iRun,1) = sum(SimParams.totUserRateE);
            SimParams.groupSumRate.tsca(iRun,1) = 1;
            SimParams.groupSumRate.isSucceded(iRun,1) = 1;
            
            gprintf('Completed [%s] with SR - %g \n',SimParams.innerPrecoder,SimParams.groupSumRate.srate(iRun,1));
        end
        
    case 'CVX'
        
        if length(precType) ~= 2
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        nSCAIterations = SimParams.SCA;
        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            alphaPowerReRun = true;
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            while alphaPowerReRun
                
                cvxInnerMP = complex(randn(SimParams.tStatBeams,SimParams.nUsers),...
                    randn(SimParams.tStatBeams,SimParams.nUsers)) / sqrt(2);
                
                cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
                    + alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'CVX');
                
%                 SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
%                 userBetaP = SimParams.totUserBeta;

                userBetaP = ones(SimParams.nUsers,1);
                
                for iSca = 1:nSCAIterations
                    
                    cvxInnerM = sdpvar(SimParams.tStatBeams,SimParams.nUsers,'full','complex');
                    userRate = sdpvar(SimParams.nUsers,1);
                    userBeta = sdpvar(SimParams.nUsers,1);
                    
                    gConstraints = [];
                    for iGroup = 1:SimParams.nGroups
                        
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM;
                            totalIF = sqrt(SimParams.N0);
                            
                            for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                                if jUser ~= iUser
                                    xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                    totalIF = [totalIF, effChannel * cvxInnerM(:,xUserIndex)];
                                end
                            end
                            
                            switch precType{2}
                                case 'Opt'
                                    for jGroup = 1:SimParams.nGroups
                                        if jGroup ~= iGroup
                                            xUserIndices = SimParams.groupInfo(jGroup).gUserIndices;
                                            totalIF = [totalIF, effChannel * cvxInnerM(:,xUserIndices)];
                                        end
                                    end
                                case 'Ignore'
                                    % nothing performed (neighboring interference is not considered
                                otherwise
                                    totalIF = [totalIF, sqrt(str2double(precType{2})) * ones(1,SimParams.nGroups - 1)];
                            end
                            
                            gConstraints = [gConstraints, rcone(totalIF(:),userBeta(cUserIndex,1),1/2)];
                            
                            pUser = real(effChannel * cvxInnerM(:,cUserIndex));
                            qUser = imag(effChannel * cvxInnerM(:,cUserIndex));
                            pBarUser = real(effChannel * cvxInnerMP(:,cUserIndex));
                            qBarUser = imag(effChannel * cvxInnerMP(:,cUserIndex));
                            
                            xSignal = (pBarUser^2 + qBarUser^2) / userBetaP(cUserIndex,1) ...
                                            + 2 * pBarUser * (pUser - pBarUser) / userBetaP(cUserIndex,1) ...
                                            + 2 * qBarUser * (qUser - qBarUser) / userBetaP(cUserIndex,1) ...
                                            - ((pBarUser^2 + qBarUser^2) / userBetaP(cUserIndex,1)^2) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                                        
                            gConstraints = [gConstraints, userRate(cUserIndex,1) == 1 + xSignal];                                                        
                        end
                        
                        if ~any(strcmpi(precType{2},{'Opt','Ignore'}))
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndices = SimParams.groupInfo(iGroup).gUserIndices;
                                        neighborIF = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * sBeamM * cvxInnerM(:,xUserIndices);
                                        gConstraints = [gConstraints, rcone(neighborIF,str2double(precType{2}),1/2)];
                                    end
                                end
                            end
                        end
                    end
                    
                    vectorizedBeams = sBeamM * cvxInnerM;
                    gConstraints = [gConstraints, rcone(vectorizedBeams(:),SimParams.txPower,1/2)];

                    optObjective = -geomean(userRate);
                    probOptions = sdpsettings('verbose',0,'solver','sedumi');
                    probSolution = optimize(gConstraints,optObjective,probOptions);
                    
                    if probSolution.problem >= 0                 
                        alphaSca = 1;
                        SimParams.sBeamM = sBeamM;
                        userRate = value(userRate);
                        cvxInnerMP = value(cvxInnerM) * alphaSca + cvxInnerMP * (1 - alphaSca);
                        userBetaP = value(userBeta) * alphaSca + userBetaP * (1 - alphaSca);
                        
                        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
                        SimParams.groupSumRate.srate(iRun,iSca) = sum(SimParams.totUserRateE);
                        SimParams.groupSumRate.tsca(iRun,1) = iSca;
                        
                        if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-1) && (iSca > 5)) || (iSca == nSCAIterations)
                            gprintf('Completed [%s] with SR - %g actual rate - [%g]\n',SimParams.innerPrecoder,log2(prod(userRate)),SimParams.groupSumRate.srate(iRun,iSca));
                            alphaPowerReRun = 0;
                            SimParams.groupSumRate.isSucceded(iRun,1) = 1;
                            break;
                        else
                            gprintf('SR progress for [%s] SCA iteration - [%3d] %f, actual rate [%g] \n',SimParams.innerPrecoder,iSca,log2(prod(userRate)),SimParams.groupSumRate.srate(iRun,iSca));
                        end
                    else
                        gprintf('Convergence failure !! \n');
                        alphaPowerReRun = 0;
                        SimParams.groupSumRate.isSucceded(iRun,1) = 0;
                        break;
                    end
                    
                end
                
            end
            
        end
        
    case 'CVXG' % Without constraining interference to neighboring groups (fixed interference values can be considered)
        
        if length(precType) ~= 2
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        nSCAIterations = SimParams.SCA;
        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            alphaPowerReRun = 1;
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            while alphaPowerReRun
                
                cvxInnerMP = initializeSCAPoints(SimParams,'CVXG');
                randoInnerMP = complex(randn(size(cvxInnerMP)),randn(size(cvxInnerMP))) / sqrt(2);
                cvxInnerMP = alphaM^alphaPowerReRun * cvxInnerMP + (1 - alphaM^alphaPowerReRun) * randoInnerMP;
                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVXG');
                userBetaP = SimParams.totUserBeta;
                
                for iSca = 1:nSCAIterations
                    
                    cvx_begin
                    
                    variable cvxInnerM(SimParams.gStatBeams,SimParams.nUsers) complex
                    variables userRate(SimParams.nUsers,1) userBeta(SimParams.nUsers,1) groupBeta(SimParams.nUsers,SimParams.nGroups)
                    expressions xSignal totalIF
                    
                    maximize geo_mean(userRate)
                    
                    subject to
                    
                    for iGroup = 1:SimParams.nGroups
                        
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                            totalIF = sqrt(SimParams.N0);
                            
                            %Interference from my own group: Intragroup
                            for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                                if jUser ~= iUser
                                    xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                    totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,xUserIndex)];
                                end
                            end
                            
                            switch precType{2}
                                case 'Opt'
                                    for jGroup = 1:SimParams.nGroups
                                        if jGroup ~= iGroup
                                            xUserIndices = SimParams.groupInfo(jGroup).gUserIndices;
                                            totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(jGroup).statBeams * cvxInnerM(:,xUserIndices)];
                                        end
                                    end
                                case 'Ignore'
                                    % nothing performed (neighboring interference is not considered
                                otherwise
                                    totalIF = [totalIF, sqrt(str2double(precType{2})) * ones(1,SimParams.nGroups - 1)];
                            end
                            
                            {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            
                            if strcmpi(precType{2},'Opt') % Assuming other groups contribute zero interference (using outer beamformer design)
                                {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            else
                                totalIF = [totalIF, sqrt(str2double(precType{2}))];
                                {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            end
                            
                            xSignal = abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1) ...
                                + 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex)) / userBetaP(cUserIndex,1) ...
                                - (abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / (userBetaP(cUserIndex,1).^2)) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                            userRate(cUserIndex,1) == 1 + real(xSignal);
                            0 == imag(xSignal);
                            
                        end
                        
                        if ~any(strcmpi(precType{2},{'Opt','Ignore'}))
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndices = SimParams.groupInfo(iGroup).gUserIndices;
                                        neighborIF = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,xUserIndices);
                                        {neighborIF.',str2double(precType{2}),1} <In> rotated_complex_lorentz(length(neighborIF));
                                    end
                                end
                            end
                        end
                        
                    end
                    
                    vectorizedBeams = [];
                    for iGroup = 1:SimParams.nGroups
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            tempBeam = SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,cUserIndex);
                            vectorizedBeams = [vectorizedBeams; tempBeam(:)];
                        end
                    end
                    
                    {vectorizedBeams,SimParams.txPower,1} <In> rotated_complex_lorentz(length(vectorizedBeams));
                    
                    cvx_end
                    
                    if contains(cvx_status,'Solved') && ~any(isnan(cvxInnerM(:)))
                        SimParams.sBeamM = sBeamM;
                        cvxInnerMP = cvxInnerM;
                        userBetaP = userBeta;
                        
                        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerM,'CVXG');
                        SimParams.groupSumRate.srate(iRun,iSca) = sum(SimParams.totUserRateE);
                        SimParams.groupSumRate.tsca(iRun,1) = iSca;
                        
                        if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-1) && (iSca > 5)) || (iSca == nSCAIterations)
                            gprintf('Completed [%s] with SR - %g actual rate - [%g]\n',SimParams.innerPrecoder,log2(prod(userRate)),SimParams.groupSumRate.srate(iRun,iSca));
                            alphaPowerReRun = 0;
                            SimParams.groupSumRate.isSucceded(iRun,1) = 1;
                            break;
                        else
                            gprintf('SR progress for [%s] SCA iteration - [%3d] %f, actual rate [%g] \n',SimParams.innerPrecoder,iSca,log2(prod(userRate)),SimParams.groupSumRate.srate(iRun,iSca));
                        end
                    else
                        gprintf('Convergence failure !! \n');
                        alphaPowerReRun = 0;
                        SimParams.groupSumRate.isSucceded(iRun,1) = 0;
                        break;
                    end
                    
                end
                
            end
            
        end
        
end

end

function nprintf(varargin)
end
