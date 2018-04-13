
function SimParams = digitalBeamformerDesign(SimParams)

if ispc
    gprintf = @fprintf;
else
    gprintf = @nprintf;
end

% cvx_quiet('true');
% cvx_expert('true');
% cvx_solver('sdpt3');

alphaM = 0.75;
sBeamM = SimParams.sBeamM;
SimParams.tStatBeams = size(sBeamM,2);
precType = strsplit(SimParams.innerPrecoder,'_');

switch precType{1}
    
    case 'ZF'
        
        for iRun = 1:SimParams.nMontRuns
            
            gprintf('Montecarlo run : [%d] with Statistical Beams per group : [%d / %d] for {%s} \n',iRun,SimParams.gStatBeams,SimParams.tStatBeams,SimParams.statBeamType);
            
            H = zeros(SimParams.nUsers,SimParams.nTransmit);
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                    H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser);
                end
            end
            
            tStart = tic;
            effChannel = H * sBeamM; % * zeroMatrix;
            ZF = effChannel' / (effChannel * effChannel');
            
            zfPower = diag(1./sqrt(diag(ZF' * ZF)));
            wfPower = SimParams.txPower * 0.9;
            %             wfPower = waterfill(SimParams.txPower,reshape(diag(ZF' * ZF),1,[]));
            ZF = ZF * sqrt(diag(wfPower)) * zfPower;
            eZF = sBeamM * ZF;
            cvxInnerMP = ZF * sqrt(SimParams.txPower) / sqrt(trace(eZF * eZF'));
            SimParams.sBeamM = sBeamM;
            %SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
            
            SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'ZF');
            
            SimParams.groupSumRate.srate(iRun,1) = sum(SimParams.totUserRateE);
            SimParams.groupSumRate.tsca(iRun,1) = 1;
            SimParams.groupSumRate.stime(iRun,1) = toc(tStart);
            SimParams.groupSumRate.isSucceded(iRun,1) = 1;
            
            gprintf('Completed [%s] with SR - %g \n',SimParams.innerPrecoder,SimParams.groupSumRate.srate(iRun,1));
        end
        
    case 'CVX'
        
        if length(precType) == 1
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        nSCAIterations = SimParams.SCA;
        SimParams.groupSumRate.sbeta = cell(SimParams.nMontRuns,1);
        
        for iRun = 1:SimParams.nMontRuns
            
            gprintf('Montecarlo run : [%d] with Statistical Beams per group : [%d / %d] for {%s} \n',iRun,SimParams.gStatBeams,SimParams.tStatBeams,SimParams.statBeamType);
            
            alphaPowerReRun = true;
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            while alphaPowerReRun
                
                tStart = tic;
                cvxInnerMP = complex(randn(SimParams.tStatBeams,SimParams.nUsers),...
                    randn(SimParams.tStatBeams,SimParams.nUsers)) / sqrt(2);
                
                cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
                    + alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'CVX');
                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
                userBetaP = SimParams.totUserBeta;
                
                for iSca = 1:nSCAIterations
                    
                    cvx_begin
                    
                    variable cvxInnerM(SimParams.tStatBeams,SimParams.nUsers) complex
                    variables userRate(SimParams.nUsers,1) userBeta(SimParams.nUsers,1)
                    expressions xSignal totalIF pUser(SimParams.nUsers,1) qUser(SimParams.nUsers,1) pBarUser(SimParams.nUsers,1) qBarUser(SimParams.nUsers,1)
                    
                    maximize geo_mean(userRate)
                    
                    subject to
                    
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
                                case 'Discard'
                                case 'Ignore'
                                    % nothing performed (neighboring interference is not considered
                                otherwise
                                    totalIF = [totalIF, sqrt(str2double(precType{2})) * ones(1,SimParams.nGroups - 1)];
                            end
                            
                            {totalIF(:),userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            pUser(cUserIndex,1) = real(effChannel * cvxInnerM(:,cUserIndex));
                            qUser(cUserIndex,1) = imag(effChannel * cvxInnerM(:,cUserIndex));
                            pBarUser(cUserIndex,1) = real(effChannel * cvxInnerMP(:,cUserIndex));
                            qBarUser(cUserIndex,1) = imag(effChannel * cvxInnerMP(:,cUserIndex));
                            
                            xSignal = (pBarUser(cUserIndex,1)^2 + qBarUser(cUserIndex,1)^2) / userBetaP(cUserIndex,1) ...
                                + 2 * pBarUser(cUserIndex,1) * (pUser(cUserIndex,1) - pBarUser(cUserIndex,1)) / userBetaP(cUserIndex,1) ...
                                + 2 * qBarUser(cUserIndex,1) * (qUser(cUserIndex,1) - qBarUser(cUserIndex,1)) / userBetaP(cUserIndex,1) ...
                                - ((pBarUser(cUserIndex,1)^2 + qBarUser(cUserIndex,1)^2) / userBetaP(cUserIndex,1)^2) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                            userRate(cUserIndex,1) == 1 + xSignal;
                        end
                        
                        if ~any(strcmpi(precType{2},{'Opt','Ignore'}))
                            if strcmpi(precType{2},'Discard')
                                ifLimit = str2double(precType{3});
                            else
                                ifLimit = str2double(precType{2});
                            end
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndices = SimParams.groupInfo(iGroup).gUserIndices;
                                        neighborIF = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * sBeamM * cvxInnerM(:,xUserIndices);
                                        {neighborIF.',ifLimit,1} <In> rotated_complex_lorentz(length(neighborIF));
                                    end
                                end
                            end
                        end
                    end
                    
                    vectorizedBeams = sBeamM * cvxInnerM;
                    {vectorizedBeams(:),SimParams.txPower,1} <In> rotated_complex_lorentz(length(vectorizedBeams(:)));
                    
                    cvx_end
                    
                    if contains(cvx_status,'Solved') && ~any(isnan(cvxInnerM(:)))
                        alphaSca = 1;
                        SimParams.sBeamM = sBeamM;
                        cvxInnerMP = cvxInnerM * alphaSca + cvxInnerMP * (1 - alphaSca);
                        userBetaP = userBeta * alphaSca + userBetaP * (1 - alphaSca);
                        
                        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerM,'CVX');
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
                
                SimParams.groupSumRate.stime(iRun,1) = toc(tStart);
                SimParams.groupSumRate.sbeta{iRun,1} = zeros(SimParams.nGroups,SimParams.nUsers);
                
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM;
                        effGroupIF = zeros(SimParams.nGroups,1);
                        for jGroup = 1:SimParams.nGroups
                            if jGroup ~= iGroup
                                xUserIndices = SimParams.groupInfo(jGroup).gUserIndices;
                                effGroupIF(jGroup,1) = norm(effChannel * cvxInnerM(:,xUserIndices))^2;
                            end
                        end
                        SimParams.groupSumRate.sbeta{iRun,1}(:,cUserIndex) = effGroupIF;
                    end
                end
                
            end
            
        end
        
    case 'CVXG'
        
        if length(precType) == 1
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        nSCAIterations = SimParams.SCA;
        for iRun = 1:SimParams.nMontRuns
            
            gprintf('Montecarlo run : [%d] with Statistical Beams per group : [%d / %d] for {%s} \n',iRun,SimParams.gStatBeams,SimParams.tStatBeams,SimParams.statBeamType);
            
            alphaPowerReRun = true;
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            while alphaPowerReRun
                
                tStart = tic;
                cvxInnerMP = complex(randn(SimParams.gStatBeams,SimParams.nUsers),...
                    randn(SimParams.gStatBeams,SimParams.nUsers)) / sqrt(2);
                
                cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
                    + alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'CVXG');
                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVXG');
                userBetaP = SimParams.totUserBeta;
                
                for iSca = 1:nSCAIterations
                    
                    cvx_begin
                    
                    variable cvxInnerM(SimParams.gStatBeams,SimParams.nUsers) complex
                    variables userRate(SimParams.nUsers,1) userBeta(SimParams.nUsers,1) groupBeta(SimParams.nUsers,SimParams.nGroups)
                    expressions xSignal totalIF pUser(SimParams.nUsers,1) qUser(SimParams.nUsers,1) pBarUser(SimParams.nUsers,1) qBarUser(SimParams.nUsers,1)
                    
                    maximize geo_mean(userRate)
                    
                    subject to
                    
                    for iGroup = 1:SimParams.nGroups
                        
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                            totalIF = sqrt(SimParams.N0);
                            
                            for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                                if jUser ~= iUser
                                    xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                    totalIF = [totalIF, effChannel * cvxInnerM(:,xUserIndex)];
                                end
                            end
                            
                            switch precType{2}
                                case 'Opt'
                                    totalIF = [totalIF, groupBeta(cUserIndex,:)];
                                case 'Ignore'
                                    % nothing performed (neighboring interference is not considered
                                otherwise
                                    totalIF = [totalIF, sqrt(str2double(precType{2})) * ones(1,SimParams.nGroups - 1)];
                            end
                            
                            {totalIF(:),userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            pUser(cUserIndex,1) = real(effChannel * cvxInnerM(:,cUserIndex));
                            qUser(cUserIndex,1) = imag(effChannel * cvxInnerM(:,cUserIndex));
                            pBarUser(cUserIndex,1) = real(effChannel * cvxInnerMP(:,cUserIndex));
                            qBarUser(cUserIndex,1) = imag(effChannel * cvxInnerMP(:,cUserIndex));
                            
                            xSignal = (pBarUser(cUserIndex,1)^2 + qBarUser(cUserIndex,1)^2) / userBetaP(cUserIndex,1) ...
                                + 2 * pBarUser(cUserIndex,1) * (pUser(cUserIndex,1) - pBarUser(cUserIndex,1)) / userBetaP(cUserIndex,1) ...
                                + 2 * qBarUser(cUserIndex,1) * (qUser(cUserIndex,1) - qBarUser(cUserIndex,1)) / userBetaP(cUserIndex,1) ...
                                - ((pBarUser(cUserIndex,1)^2 + qBarUser(cUserIndex,1)^2) / userBetaP(cUserIndex,1)^2) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                            userRate(cUserIndex,1) == 1 + xSignal;
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
                        
                        if any(strcmpi(precType{2},{'Opt'}))
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        jUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        neighborIF = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,SimParams.groupInfo(iGroup).gUserIndices);
                                        {neighborIF.',groupBeta(jUserIndex,iGroup),1} <In> rotated_complex_lorentz(length(neighborIF));
                                    end
                                end
                            end
                        end
                        
                    end
                    
                    for iGroup = 1:SimParams.nGroups
                        vectorizedBeams = SimParams.groupInfo(iGroup).statBeams * cvxInnerM(:,SimParams.groupInfo(iGroup).gUserIndices);
                        {vectorizedBeams(:),(SimParams.txPower / SimParams.nGroups),1} <In> rotated_complex_lorentz(length(vectorizedBeams(:)));
                    end
                    
                    cvx_end
                    
                    if contains(cvx_status,'Solved') && ~any(isnan(cvxInnerM(:)))
                        alphaSca = 1;
                        SimParams.sBeamM = sBeamM;
                        cvxInnerMP = cvxInnerM * alphaSca + cvxInnerMP * (1 - alphaSca);
                        userBetaP = userBeta * alphaSca + userBetaP * (1 - alphaSca);
                        
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
                
                SimParams.groupSumRate.stime(iRun,1) = toc(tStart);
                SimParams.groupSumRate.sbeta{iRun,:} = groupBeta;
                
            end
            
        end
        
    case 'KKT'
        
        stInstant = 1;
        cKG = 5 * ones(SimParams.nUsers,SimParams.nGroups);
        if length(precType) == 3
            cKG = str2double(precType{1,3}) * ones(SimParams.nUsers,SimParams.nGroups);
        end
        
        bK = zeros(SimParams.nUsers,1);
        userGamma = zeros(SimParams.nUsers,1);
        
        if length(precType) == 1
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        %Generate the user specific channel
        for iGroup = 1:SimParams.nGroups
            SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
            end
        end
        
        cvxInnerMP = complex(randn(SimParams.gStatBeams,SimParams.nUsers),...
            randn(SimParams.gStatBeams,SimParams.nUsers)) / sqrt(2);
        
        %cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
            %+ alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'KKT');
        
        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVXG');
        userBetaP = SimParams.totUserBeta;
                
        cvxInnerM = zeros(size(cvxInnerMP));
        userBeta = zeros(size(userBetaP));
        
        %aK = ones(SimParams.nUsers,1) * 0.5;
        
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                userGammaP(cUserIndex,1) = norm( effChannel * cvxInnerMP(:,cUserIndex)).^2 / userBetaP(cUserIndex,1);
            end
        end
        
        for iIterations = stInstant:SimParams.innerIterations
            
            for interIterations = 1:SimParams.interIterations
                
                aK = log2(exp(1))./(1+userGammaP);
                
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        bK(cUserIndex,1) = aK(cUserIndex,1) * norm(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1).^2;
                    end
                end
                
                intermediateY = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nUsers);
                intermediateZ = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nGroups);
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if iUser ~= jUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                intermediateY(:,:,cUserIndex) = intermediateY(:,:,cUserIndex) + effChannel' * effChannel * bK(xUserIndex,1);
                            end
                        end
                    end
                    
                    switch precType{2}
                        case 'Cen'
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        yUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                        intermediateY(:,:,cUserIndex) = intermediateY(:,:,cUserIndex) + effChannel' * effChannel * bK(yUserIndex,1);
                                    end
                                end
                            end
                        case 'Opt'
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        yUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                        intermediateZ(:,:,iGroup) = intermediateZ(:,:,iGroup) + effChannel' * effChannel * cKG(yUserIndex,iGroup);
                                    end
                                end
                            end
                        case 'Ignore'
                            % nothing performed (neighboring interference is not considered
                        otherwise
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                        intermediateZ(:,:,iGroup) = intermediateZ(:,:,iGroup) + str2double(precType{1,2}) * eye(SimParams.gstatBeams) * cKG(xUserIndex,iGroup);
                                    end
                                end
                            end
                    end
                end
                
                for iGroup = 1:SimParams.nGroups
                    
                    wTemp = zeros(SimParams.groupInfo(iGroup).nUsers,SimParams.gStatBeams);
                    wOverall = zeros(SimParams.groupInfo(iGroup).nUsers,SimParams.gStatBeams);
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        wTemp(iUser,:) = (aK(cUserIndex,1) * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) ) / userBetaP(cUserIndex,1);
                    end
                    
                    muMax = 1e5;
                    muMin = 0;
                    iterateAgain = 1;
                    
                    while iterateAgain
                        currentMu = (muMax + muMin) / 2;
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            wOverall(iUser,:) = wTemp(iUser,:) * pinv(intermediateY(:,:,cUserIndex) + intermediateZ(:,:,iGroup) + (SimParams.groupInfo(iGroup).statBeams' * SimParams.groupInfo(iGroup).statBeams) * currentMu);
                        end
                        
                        totalPower = real(trace(SimParams.groupInfo(iGroup).statBeams * (wOverall' * wOverall) * SimParams.groupInfo(iGroup).statBeams'));
                        
                        if totalPower > (SimParams.txPower / SimParams.nGroups)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if (muMax - muMin) <= 1e-3
                            iterateAgain = 0;
                        end
                    end
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        cvxInnerM(:,cUserIndex) = wOverall(iUser,:)';
                    end
                    
                end
                
                
                %%%%%%%%%%%% Calculate Interference term beta
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        userBeta(cUserIndex,1) = SimParams.N0;
                        
                        %Intra-Group Interference calculation
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if jUser ~= iUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                userBeta(cUserIndex,1) = userBeta(cUserIndex,1) + norm(effChannel * cvxInnerM(:,xUserIndex))^2;
                            end
                        end
                        
                        %Inter-Group Interference Calculation
                        for jGroup = 1:SimParams.nGroups
                            if jGroup ~= iGroup
                                switch precType{2}
                                    case 'Opt'
                                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(jGroup).statBeams;
                                        userBeta(cUserIndex,1) = userBeta(cUserIndex,1) + norm(effChannel * cvxInnerM(:,SimParams.groupInfo(jGroup).gUserIndices))^2;
                                    case 'Ignore'
                                        % nothing performed (neighboring interference is not considered
                                    otherwise
                                        userBeta(cUserIndex,1) = userBeta(cUserIndex,1) +  str2double(precType{2});
                                end
                            end
                        end
                    end
                    
                end
                
                %Update Gamma
                
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        userGamma(cUserIndex,1) = norm( effChannel * cvxInnerM(:,cUserIndex)).^2 / userBeta(cUserIndex,1);
                    end
                end
                
                % Updating for next iteration
                
                userBetaP = userBeta;
                cvxInnerMP = cvxInnerM;
                userGammaP = real(userGamma) .* real(userGamma>0);

                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerM,'CVXG');
                SimParams.groupSumRate.srate(iIterations,interIterations) = sum(SimParams.totUserRateE);
                SimParams.groupSumRate.tsca(iIterations,1) = interIterations;
                SimParams.groupSumRate.isSucceded(iIterations,1) = 1;
                
            end
        end
        
        
        %end
        %end
        
        
    case 'KKTU'
        
        stInstant = 1;
        cKG = 10 * ones(SimParams.nUsers,SimParams.nGroups);
        
        aK = zeros(SimParams.nUsers,1);
        bK = zeros(SimParams.nUsers,1);
        
        if length(precType) == 1
            precType{2} = 'Opt'; % No interference is considered from neighboring group (0) else value can be specified (1 - Inf) (fixed power)
        end
        
        alphaPowerReRun = true;
        %Generate the user specific channel
        for iGroup = 1:SimParams.nGroups
            SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
            end
        end
        
        cvxInnerMP = complex(randn(SimParams.gStatBeams,SimParams.nUsers),...
            randn(SimParams.gStatBeams,SimParams.nUsers)) / sqrt(2);
        
        %cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
           % + alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'CVXG');
        
        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVXG');
        userBetaP = SimParams.totUserBeta;
        %userGammaP = SimParams.totuserGamma;
        
        cvxInnerM = zeros(size(cvxInnerMP));
        userBeta = zeros(size(userBetaP));
        
        for iGroup = 1:SimParams.nGroups
            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                userGammaP(cUserIndex,1) = norm( effChannel * cvxInnerMP(:,cUserIndex)).^2 / userBetaP(cUserIndex,1);
            end
        end
        
        userGamma = zeros(size(userGammaP));
        %aK = ones(SimParams.nUsers,1) * 0.9;
        
        %aK = log2(exp(1))./(1+userGammaP);
        
        for iIterations = stInstant:SimParams.innerIterations
            
            phiK = sqrt(userGammaP ./ userBetaP);
            
            for interIterations = 1:SimParams.interIterations
                
                
                
                aK = log2(exp(1)) * 2 * phiK ./ (1 + userGammaP);
                
                %                 for iGroup = 1:SimParams.nGroups
                %                     for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                %                         cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                %                         effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                %                         bK(cUserIndex,1) = aK(cUserIndex,1) * norm(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1).^2;
                %                     end
                %                 end
                
                %bK = max(bK,0);
                
                bK = 0.5 * aK .* phiK;
                
                intermediateY = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nUsers);
                intermediateZ = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nGroups);
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if iUser ~= jUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                intermediateY(:,:,cUserIndex) = intermediateY(:,:,cUserIndex) + effChannel' * effChannel * bK(xUserIndex,1);
                            end
                        end
                    end
                    
                    switch precType{2}
                        case 'Opt'
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        yUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                        intermediateY(:,:,iGroup) = intermediateY(:,:,iGroup) + effChannel' * effChannel * cKG(yUserIndex,iGroup);
                                    end
                                end
                            end
                        case 'Ignore'
                            % nothing performed (neighboring interference is not considered
                        otherwise
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                        intermediateZ(:,:,iGroup) = intermediateZ(:,:,iGroup) + effChannel' * effChannel * cKG(xUserIndex,iGroup);
                                    end
                                end
                            end
                    end
                end
                
                for iGroup = 1:SimParams.nGroups
                    
                    wTemp = zeros(SimParams.groupInfo(iGroup).nUsers,SimParams.gStatBeams);
                    wOverall = zeros(SimParams.groupInfo(iGroup).nUsers,SimParams.gStatBeams);
                    
%                     for iUser = 1:SimParams.groupInfo(iGroup).nUsers
%                         cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
%                         effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
%                         wTemp(iUser,:) = (aK(cUserIndex,1) * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) )./ userBetaP(cUserIndex,1);
%                     end
                    
                    muMax = 100000;
                    muMin = 0;
                    iterateAgain = 1;
                    
                    while iterateAgain
                        currentMu = (muMax + muMin) / 2;
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                            wOverall(iUser,:) = 0.5 *  aK(cUserIndex,1) * (pinv(intermediateY(:,:,cUserIndex) + intermediateZ(:,:,iGroup) + (SimParams.groupInfo(iGroup).statBeams' * SimParams.groupInfo(iGroup).statBeams) * currentMu)) * effChannel';
                            %wOverall(iUser,:) =  pinv(intermediateY(:,:,cUserIndex) + intermediateZ(:,:,iGroup) + (SimParams.groupInfo(iGroup).statBeams' * SimParams.groupInfo(iGroup).statBeams) * currentMu) * wTemp(iUser,:)';
                        end
                        
                        totalPower = real(trace(SimParams.groupInfo(iGroup).statBeams * (wOverall' * wOverall) * SimParams.groupInfo(iGroup).statBeams'));
                        
                        if totalPower > (SimParams.txPower / SimParams.nGroups)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if (muMax - muMin) <= 1e-6
                            iterateAgain = 0;
                        end
                    end
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        cvxInnerM(:,cUserIndex) = wOverall(iUser,:)';
                    end
                    
                end
                
                
                %%%%%%%%%%%% Calculate Interference term beta
                
                for iGroup = 1:SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        userBeta(cUserIndex,1) = SimParams.N0;
                        %Intra-Group Interference calculation
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if jUser ~= iUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                userBeta(cUserIndex,1) = userBeta(cUserIndex,1) + norm(effChannel * cvxInnerM(:,xUserIndex))^2;
                            end
                        end
                        
                        for jGroup = 1:SimParams.nGroups
                            if jGroup ~= iGroup
                                switch precType{2}
                                    case 'Opt'
                                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(jGroup).statBeams;
                                        userBeta(cUserIndex,1) = userBeta(cUserIndex,1) + norm(effChannel * cvxInnerM(:,SimParams.groupInfo(jGroup).gUserIndices))^2;
                                    case 'Ignore'
                                        % nothing performed (neighboring interference is not considered
                                    otherwise
                                        userBeta(cUserIndex,1) = userBeta(cUserIndex,1) +  str2double(precType{2});
                                end
                            end
                        end
                    end
                    
                end
                
                %Update Gamma
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        userGamma(cUserIndex,1) = (abs(effChannel * cvxInnerMP(:,cUserIndex)) - 0.5 * phiK(cUserIndex,1) * userBetaP(cUserIndex,1)) * 2 ./ phiK(cUserIndex,1);
                        %userGamma(cUserIndex,1)  = userGamma(cUserIndex,1)  + prevGamma * 0.1;
                        %prevGamma = 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex))  ./ userBetaP(cUserIndex,1);
                        %userGamma(cUserIndex,1) = prevGamma + (norm(effChannel * cvxInnerMP(:,cUserIndex))^2 ./ userBetaP(cUserIndex,1)) * (1 - ((userBeta(cUserIndex,1) - userBetaP(cUserIndex,1)) / userBetaP(cUserIndex,1)));
                    end
                end
                
                % Updating for next iteration

                userGammaP = real(userGamma) .* (real(userGamma) > 0);

                userBetaP = userBeta;
                cvxInnerMP = cvxInnerM;
                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerM,'CVXG');
                SimParams.groupSumRate.srate(iIterations,interIterations) = sum(SimParams.totUserRateE);
                SimParams.groupSumRate.tsca(iIterations,1) = interIterations;
                SimParams.groupSumRate.isSucceded(iIterations,1) = 1;
                
            end
        end
end

%end
%end




end



function nprintf(varargin)
end
