
function SimParams = digitalBeamformerDesign(SimParams)

if ispc
    gprintf = @fprintf;
else
    gprintf = @nprintf;
end

cvx_quiet('true');
cvx_expert('true');
cvx_solver('sdpt3');

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
        stepFactor = 1;
        stInstant = 1;
        incCounter = 1;
        
        cK = 1e-1 * ones(SimParams.nUsers,1);
        
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
        
        %while alphaPowerReRun
        
        %tStart = tic;
        cvxInnerMP = complex(randn(SimParams.gStatBeams,SimParams.nUsers),...
            randn(SimParams.gStatBeams,SimParams.nUsers)) / sqrt(2);
        
        cvxInnerMP = (1 - alphaM^alphaPowerReRun) * cvxInnerMP ...
            + alphaM^alphaPowerReRun * initializeSCAPoints(SimParams,'KKT');
        
        SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'KKT');
        userBetaP = SimParams.totUserBeta;
        userGammaP = SimParams.totuserGamma;
        
        for iIterations = stInstant:SimParams.innerIterations
            
            stepG = stepFactor;
            
            for interIterations = 1:SimParams.interIterations
                
                stepG = stepG * 0.9;
                aK = 1./(1+userGammaP); %Dual variable aK
                
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        realRatio(cUserIndex,1) = real(effChannel * cvxInnerMP(:,cUserIndex));
                        sqVal(cUserIndex,1) = realRatio(cUserIndex,1).^2;
                    end
                end
                intermediateX = sqVal./(userBetaP).^2;
                bK = aK .* intermediateX; %Dual variable bK
                
                intermediateY = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nUsers);
                intermediateZ = zeros(SimParams.gStatBeams,SimParams.gStatBeams,SimParams.nUsers);
                
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        intermediateY(:,:,iUser) = zeros(SimParams.gStatBeams,SimParams.gStatBeams);
                        xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChan = SimParams.groupInfo(iGroup).userChannel(:,:,xUserIndex) * SimParams.groupInfo(iGroup).statBeams;
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if iUser ~= jUser
                                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;
                                intermediateY(:,:,iUser) = intermediateY(:,:,iUser) + effChannel'*effChannel * bK(jUser,1);
                            end
                        end
                    end
                    
                    if any(strcmpi(precType{2},{'Opt'}))
                        for jGroup = 1:SimParams.nGroups
                            if jGroup ~= iGroup
                                for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                    %jUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                    effChannel = SimParams.groupInfo(jGroup).userChannel(:,:,jUser) * SimParams.groupInfo(iGroup).statBeams;% * cvxInnerMP(:,SimParams.groupInfo(iGroup).gUserIndices);
                                    intermediateZ(:,:,jUser) = intermediateZ(:,:,jUser) + effChannel'*effChannel * cK(jUser,1);
                                    %{neighborIF.',groupBeta(jUserIndex,iGroup),1} <In> rotated_complex_lorentz(length(neighborIF));
                                end
                            end
                        end
                    end
                end
                
                for iGroup = 1:SimParams.nGroups
                    muMax = 100000;
                    muMin = 0;
                    iterateAgain = 1;
                    nBands = 1;
                    while iterateAgain
                        totalPower = 0;
                        currentMu = (muMax + muMin) / 2;
                        for iBand = 1:nBands
                            for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                                effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                                cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                                
                                %Dual variable update w
                                %interVal(:,cUserIndex) = aK(cUserIndex,1) * (effChannel*effChannel' * cvxInnerMP(:,cUserIndex))./userBetaP(cUserIndex,1);
                                cvxInnerM(:,cUserIndex) = aK(cUserIndex,1) * ( pinv(currentMu * eye(SimParams.nTransmit) + intermediateY(:,:,cUserIndex)) + intermediateZ(:,:,cUserIndex)) * (effChannel*effChannel' * cvxInnerMP(:,cUserIndex)) ./userBetaP(cUserIndex,1);
                                totalPower = totalPower + real(trace(cvxInnerM(:,cUserIndex) * cvxInnerM(:,cUserIndex)'));
                            end
                        end
                        
                        if totalPower > sum(SimParams.txPower)
                            muMin = currentMu;
                        else
                            muMax = currentMu;
                        end
                        
                        if abs(SimParams.txPower / SimParams.nGroups - sum(SimParams.txPower)) <= 1e-6
                            iterateAgain = 0;
                        end
                    end
                    
                end
                
                
                
                %%%%%%%%%%%% Calculate Interference
                for iGroup = SimParams.nGroups
                    
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        totalIF = sqrt(SimParams.N0);
                        for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                            if jUser ~= iUser
                                xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                totalIF = [totalIF, abs(effChannel * cvxInnerM(:,xUserIndex))^2];
                            end
                        end
                        
                        switch precType{2}
                            case 'Opt'
                                totalIF = [totalIF, userBetaP(cUserIndex,:)];
                            case 'Ignore'
                                % nothing performed (neighboring interference is not considered
                            otherwise
                                totalIF = [totalIF, sqrt(str2double(precType{2})) * ones(1,SimParams.nGroups - 1)];
                        end
                        
                        %userBetaP = totalIF;
                        
                    end
                end
                %Update Gamma
                for iGroup = 1:SimParams.nGroups
                    for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                        cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                        effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                        prevGamma = (2 * (effChannel'* effChannel * cvxInnerMP(:,cUserIndex))) .* (cvxInnerM - cvxInnerMP)  ./ userBetaP(cUserIndex,1);
                        userGammaP(iUser,1) = prevGamma(iUser,1) + ((norm(effChannel * cvxInnerMP)^2 ./ userBetaP(cUserIndex,1)) .* (1 - ( (totalIF(1,cUserIndex)' - userBetaP(cUserIndex,1))./userBetaP(cUserIndex,1) )  ));
                    end
                end
                
                userGammaP = real(userGammaP) .* (real(userGammaP)>0);
                targetRate(incCounter,1) = sum(log(1+userGammaP));
                userRate(incCounter,1) = log(1+userGammaP.');
                incCounter = incCounter + 1;
                SimParams.groupSumRate.srate(incCounter,1) = targetRate(incCounter,1);
                
            end
            
        end
end

%end
%end




end



function nprintf(varargin)
end
