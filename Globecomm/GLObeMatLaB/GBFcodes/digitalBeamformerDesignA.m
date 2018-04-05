
function SimParams = digitalBeamformerDesignA(SimParams)

cvx_clear;
cvx_quiet('true');
cvx_expert('true');

maxAttempts = 5;
sBeamM = SimParams.sBeamM;
precType = strsplit(SimParams.innerPrecoder,'_');

SimParams.groupSumRate.srate = zeros(SimParams.nMontRuns,SimParams.SCA);

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
            wfPower = waterfill(SimParams.txPower,reshape(diag(ZF' * ZF),1,[]));
            ZF = ZF * sqrt(diag(wfPower)) * zfPower;
            eZF = sBeamM * ZF;
            cvxInnerMP = ZF * sqrt(SimParams.txPower) / sqrt(trace(eZF * eZF'));
            SimParams.sBeamM = sBeamM;
            SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
            
            SimParams.groupSumRate.srate(1,1) = sum(SimParams.totUserRateE);
            SimParams.groupSumRate.tsca(1,1) = 1;
            fprintf('Completed [%s] with SR - %f \n',SimParams.innerPrecoder,sum(SimParams.totUserRateE));
            
        end
        
    case 'CVX'
        
        nSCAIterations = SimParams.SCA;
        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            bRun.rate = 0;
            for xRuns = 1:maxAttempts
                
                cvxInnerMP = complex(randn(SimParams.tStatBeams,SimParams.nUsers),...
                    randn(SimParams.tStatBeams,SimParams.nUsers)) / sqrt(2);
                
                alphaM = 0.5;
                cvxInnerMP = (1 - alphaM) * cvxInnerMP + alphaM * initializeSCAPoints(SimParams,'CVX');
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVX');
                userBetaP = SimParams.totUserBeta;                                
                
                for iSca = 1:nSCAIterations
                    
                    cvx_begin
                    
                    variable cvxInnerM(SimParams.tStatBeams,SimParams.nUsers) complex
                    variables userRate(SimParams.nUsers,1) userGamma(SimParams.nUsers,1) userBeta(SimParams.nUsers,1)
                    expressions xSignal totalIF
                    
                    maximize(geo_mean(userRate))
                    
                    subject to
                    
                    for iGroup = 1:SimParams.nGroups
                        
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            totalIF = sqrt(SimParams.N0);
                            
                            for jUser = 1:SimParams.groupInfo(iGroup).nUsers
                                if jUser ~= iUser
                                    xUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,jUser);
                                    totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex)];
                                end
                            end
                            
                            for jGroup = 1:SimParams.nGroups
                                if jGroup ~= iGroup
                                    for jUser = 1:SimParams.groupInfo(jGroup).nUsers
                                        xUserIndex = SimParams.groupInfo(jGroup).gUserIndices(1,jUser);
                                        totalIF = [totalIF, SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM * cvxInnerM(:,xUserIndex)];
                                    end
                                end
                            end
                            
                            % norm(totalIF,2) <= sqrt(userBeta(cUserIndex,1));
                            {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            
                            userRate(cUserIndex,1) == 1 + userGamma(cUserIndex,1);
                            
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * sBeamM;
                            xSignal = abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1) ...
                                + 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex)) / userBetaP(cUserIndex,1) ...
                                - (abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / (userBetaP(cUserIndex,1).^2)) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                            userGamma(cUserIndex,1) == real(xSignal);
                            0 == imag(xSignal);
                            
                        end
                        
                        if SimParams.limitToGroupBeamsOnly
                            tempVector = cvxInnerM(~SimParams.groupInfo(iGroup).activeBeams,SimParams.groupInfo(iGroup).gUserIndices);
                            norm(tempVector(:)) <= 0;
                        end
                        
                    end
                    
                    vectorizedBeams = [];
                    for iGroup = 1:SimParams.nGroups
                        for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                            cUserIndex = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                            tempBeam = sBeamM * cvxInnerM(:,cUserIndex);
                            vectorizedBeams = [vectorizedBeams; tempBeam(:)];
                        end
                    end
                    
                    {vectorizedBeams,SimParams.txPower,1} <In> rotated_complex_lorentz(length(vectorizedBeams));
                    
                    cvx_end
                    
                    if contains(cvx_status,'Solved')
                        SimParams.sBeamM = sBeamM;
                        cvxInnerMP = cvxInnerM;
                        userBetaP = userBeta;
                        
                        SimParams.groupSumRate.srate(iRun,iSca) = log2(prod(userRate));
                        SimParams.groupSumRate.tsca(iRun,1) = iSca;
                        
                        if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-2) && (iSca > 10))
                            SimParams = getIntraGroupIF(SimParams,sBeamM,cvxInnerMP);
                            SimParams = getInterGroupIF(SimParams,sBeamM,cvxInnerMP);
                            fprintf('Completed [%s] with SR - %f \n',SimParams.innerPrecoder,log2(prod(userRate)));
                            break;
                        else
                            fprintf('SR progress for [%s] SCA iteration - [%3d] %f \n',SimParams.innerPrecoder,iSca,log2(prod(userRate)));
                        end
                        
                        bestRandomization.cvxInnerM = cvxInnerM;
                        bestRandomization.rate = SimParams.groupSumRate.srate(iRun,iSca);
                        
                    else
                        bestRandomization.rate = nan;
                        fprintf('Continuing with next random value \n');
                        break;
                    end
                    
                end
                
                if ~isnan(bestRandomization.rate) && (bestRandomization.rate > bRun.rate)
                    bRun = bestRandomization;
                    bRun.SimParams = SimParams;
                end
                
            end
            
            SimParams.groupSumRate.srate(iRun,:) = bRun.SimParams.groupSumRate.srate(iRun,:);
            SimParams.groupSumRate.tsca(iRun,1) = bRun.SimParams.groupSumRate.tsca(iRun,1);
            
        end
        
    case 'CVXG'
        
        if length(precType) ~= 2
            precType{2} = 'Opt';
        end
        
        nSCAIterations = SimParams.SCA;
        SimParams.tStatBeams = size(sBeamM,2);
        
        for iRun = 1:SimParams.nMontRuns
            
            %Generate the user specific channel
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                end
            end
            
            bRun.rate = 0;
            for xRuns = 1:maxAttempts
                
                alphaM = 0.5;
                cvxInnerMP = initializeSCAPoints(SimParams,'CVXG');
                randoInnerMP = complex(randn(size(cvxInnerMP)),randn(size(cvxInnerMP))) / sqrt(2);
                cvxInnerMP = alphaM * cvxInnerMP + (1 - alphaM) * randoInnerMP;
                
                SimParams = evaluateUserRatesWithGroupPrecoders(SimParams,cvxInnerMP,'CVXG');
                userBetaP = SimParams.totUserBeta;
                
                for iSca = 1:nSCAIterations
                    
                    cvx_begin
                    
                    variable cvxInnerM(SimParams.gStatBeams,SimParams.nUsers) complex
                    variables userRate(SimParams.nUsers,1) userGamma(SimParams.nUsers,1) userBeta(SimParams.nUsers,1) groupBeta(SimParams.nUsers,SimParams.nGroups)
                    expressions xSignal totalIF
                    
                    maximize(geo_mean(userRate))
                    
                    subject to
                    
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
                            
                            if strcmpi(precType{2},'Opt')
                                {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            else
                                totalIF = [totalIF, sqrt(str2double(precType{2}) * (SimParams.nGroups - 1))];
                                {totalIF.',userBeta(cUserIndex,1),1} <In> rotated_complex_lorentz(length(totalIF));
                            end
                            
                            userRate(cUserIndex,1) == 1 + userGamma(cUserIndex,1);
                            effChannel = SimParams.groupInfo(iGroup).userChannel(:,:,iUser) * SimParams.groupInfo(iGroup).statBeams;
                            
                            xSignal = abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / userBetaP(cUserIndex,1) ...
                                + 2 * cvxInnerMP(:,cUserIndex)' * (effChannel' * effChannel) * (cvxInnerM(:,cUserIndex) - cvxInnerMP(:,cUserIndex)) / userBetaP(cUserIndex,1) ...
                                - (abs(effChannel * cvxInnerMP(:,cUserIndex))^2 / (userBetaP(cUserIndex,1).^2)) * (userBeta(cUserIndex,1) - userBetaP(cUserIndex,1));
                            userGamma(cUserIndex,1) == real(xSignal);
                            0 == imag(xSignal);
                            
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
                    
                    if contains(cvx_status,'Solved')
                        SimParams.sBeamM = sBeamM;
                        cvxInnerMP = cvxInnerM;
                        userBetaP = userBeta;
                        SimParams.groupSumRate.srate(iRun,iSca) = log2(prod(userRate));
                        SimParams.groupSumRate.tsca(iRun,1) = iSca;
                        
                        for iGroup = 1:SimParams.nGroups
                            SimParams.groupInfo(iGroup).Perf(iRun,:) = userRate(SimParams.groupInfo(iGroup).gUserIndices);
                            SimParams.groupInfo(iGroup).InnerBF(:,:,iRun) = cvxInnerMP(:,SimParams.groupInfo(iGroup).gUserIndices);
                        end
                        
                        if ((std(SimParams.groupSumRate.srate(iRun,max(iSca - 4,1):iSca)) <= 1e-2) && (iSca > 10))
                            fprintf('Completed [%s] with SR - %f \n',SimParams.innerPrecoder,log2(prod(userRate)));
                            break;
                        else
                            fprintf('SR progress for [%s] SCA iteration - [%3d] %f \n',SimParams.innerPrecoder,iSca,log2(prod(userRate)));
                        end
                        
                        bestRandomization.cvxInnerM = cvxInnerM;
                        bestRandomization.rate = SimParams.groupSumRate.srate(iRun,iSca);
                        
                    else
                        bestRandomization.rate = nan;
                        fprintf('Continuing with next random value \n');
                        break;
                    end
                    
                end
                
                if ~isnan(bestRandomization.rate) && (bestRandomization.rate > bRun.rate)
                    bRun = bestRandomization;
                    bRun.SimParams = SimParams;
                end
                
            end
            
            SimParams.groupSumRate.srate(iRun,:) = bRun.SimParams.groupSumRate.srate(iRun,:);
            SimParams.groupSumRate.tsca(iRun,1) = bRun.SimParams.groupSumRate.tsca(iRun,1);
            
        end
        
        
    case 'CVXSOC'
        
        K = SimParams.nUsers;
        B = sBeamM;
        S = size(sBeamM,2);
        SimParams.tStatBeams = size(sBeamM,2);
        alpha = ones(SimParams.nUsers,1);
        
        max_iter = SimParams.nRealizations;
        tolerance = 1e-3;
        
        for iSca = 1:SimParams.SCA
            
            n = 1;
            
            %Generate the user specific channel
            %H = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.nUsers);
            
            for iGroup = 1:SimParams.nGroups
                SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers);
                for iUser = 1:SimParams.groupInfo(iGroup).nUsers
                    cUser = SimParams.groupInfo(iGroup).gUserIndices(1,iUser);
                    SimParams.groupInfo(iGroup).userChannel(:,:,iUser) = getRingChannel(SimParams,iGroup,iUser,SimParams.chnType); %user specific channel
                    H(cUser,:) = SimParams.groupInfo(iGroup).userChannel(:,:,iUser); %user specific channel
                end
            end
            
            zfPrecoder = initializeSCAPoints(SimParams,'CVXSOC');
            
            signal = H * B * zfPrecoder; % Signals of all users as seen by all users
            signal_power = diag(abs(signal).^2); % Signal powers
            interference_signals = signal - diag(diag(signal));
            interference_power = sum((abs(interference_signals).^2),2); % Interference powers
            SINR = signal_power./(interference_power + 1);
            t = (1 + SINR).^alpha;
            
            % Linearization point of p,q,beta
            p_tilde = real(diag(signal));
            q_tilde = imag(diag(signal));
            beta_tilde = (p_tilde.^2 + q_tilde.^2)./(t-1);
            
            precoders(:,:,n) = zfPrecoder; % Store the first set of precoders and rate
            rate(n) = log2(prod(t));
            %                 keyboard;
            diff = 1;
            while (n <= max_iter) && (diff > tolerance)
                
                cvx_begin quiet
                
                variable w(S,K) complex
                variable t(K)
                variable beeta(K)
                expression p(K)
                expression q(K)
                expression Imat(S,K-1,K)
                
                for k=1:K % Setting up expressions
                    interferers = [1:K];
                    interferers(k) = [];
                    Imat(:,:,k) = zfPrecoder(:,interferers); % Matrices for interfering users' beamformers
                    p(k) = real(H(k,:)*B*zfPrecoder(:,k)); % Signal of user k
                    q(k) = imag(H(k,:)*B*zfPrecoder(:,k));
                end
                
                
                maximize geo_mean(t)
                
                subject to
                
                for k=1:K % Per-user interference SOC constraints
                    norm([H(k,:)*B*Imat(:,:,k) 1 0.5*(beeta(k) - 1)],2) <= 0.5*(beeta(k)+1) % SOC for interference & noise
                end
                
                for k=1:K % Linearized constraint on t
                    t(k)^(1/alpha(k)) <= 1 + 2*(p_tilde(k)/beta_tilde(k))*(p(k) - p_tilde(k)) ...
                        + 2*(q_tilde(k)/beta_tilde(k))*(q(k) - q_tilde(k)) ...
                        + ((p_tilde(k)^2 + q_tilde(k)^2)/beta_tilde(k))*(1 - ((beeta(k) - beta_tilde(k))/beta_tilde(k)))
                end
                
                norm(B*zfPrecoder, 'fro') <= sqrt(SimParams.txPower) % Total power constraint
                
                cvx_end
                rate_old = rate(n);
                n = n+1;
                rate(n) = log2(prod(t));
                %                     keyboard;
                diff = abs(rate(n) - rate_old);
                precoders(:,:,n) = zfPrecoder;
                p_tilde = p;
                q_tilde = q;
                beta_tilde = beeta;
                [SimParams.txSINR S; SimParams.nRealizations rate(n);]
                
                storage_rate(n,iSca,SimParams.nRealizations) = rate(n);
            end
            [retry_rate(iSca),ind] = max(rate); % Store the best achieved rate and precoder for the current precoder initialization
            iteration_count(iSca,SimParams.nRealizations) = length(rate);
            retry_prec(:,:,iSca) = precoders(:,:,ind);
        end
        [final_rate(SimParams.nRealizations),ind] = max(retry_rate); % Choose the maximum rate achieved among different initializations
        final_prec(:,:,SimParams.nRealizations) = retry_prec(:,:,ind);
        final_rate_zf(SimParams.nRealizations) = sumrate_zf;
        final_prec_zf(:,:,SimParams.nRealizations) = w_zf;
        final_peruserrate_zf(:,SimParams.nRealizations) = peruserrate_zf;
        final_minrate_zf(SimParams.nRealizations) = minrate_zf;
        final_userrate(:,SimParams.nRealizations) = log2(t);
        %         catch
        %             error_ch(:,:,ch) = chan;
        %             error_B(:,:,ch) = B;
        %             continue
        %         end
        for k=1:K
            userpowers(k,SimParams.nRealizations) = norm(precoders(:,k,n),2)^2;
        end
        userpowers = userpowers';
        
        
end

end