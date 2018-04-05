

function TLP(inParams,folderName)

for iBeamPerGroup = inParams.nBeamsOverGroups % parfor

    rng('default');
    SimParams.legendName = inParams.simulationLegend;
    SimParams.usersPerGroup = inParams.usersPerGroup;
    
    SimParams.nReceive = inParams.nReceive;
    SimParams.nTransmit = inParams.nTransmit;
    
    SimParams.chnType = inParams.ChannelType;
    SimParams.SCA = inParams.nSCA;
    
    SimParams.numScatterers = 20;
    SimParams.fixedGroups = 1024;
    SimParams.nRealizations = 100;
    
    SimParams.N0 = 1;
    SimParams.txSINR = inParams.txSNR;
    SimParams.txPower = 10^(SimParams.txSINR/10);
    SimParams.nUsers = sum(SimParams.usersPerGroup);

    SimParams.pathLoss = ones(1,SimParams.nUsers);
    
    SimParams.nMontRuns = inParams.nDrops;
    SimParams.statBeamType = inParams.statBeamType;
    SimParams.innerPrecoder = inParams.innerPrecoderType;
    SimParams.nGroups = length(SimParams.usersPerGroup);
    SimParams.beamsPerGroup = iBeamPerGroup / SimParams.nGroups;
    
    SimParams.totStatBeams = iBeamPerGroup;
    SimParams.gStatBeams = floor(iBeamPerGroup / SimParams.nGroups);
    SimParams.groupUserIndices = cell(SimParams.nGroups,1);
    
    SimParams.limitToGroupBeamsOnly = inParams.limitToGroupBeamsOnly;
 
    switch SimParams.chnType
        case 'Ring' %Define channel while considering ULA       
            SimParams.fixedGroupLocs = linspace(-89,90,SimParams.fixedGroups + 1);
            SimParams.fixedGroupLocs = SimParams.fixedGroupLocs(2:2:SimParams.fixedGroups);
            userClustersPerGroup = ceil(length(SimParams.fixedGroupLocs) / SimParams.nGroups);
            
        case 'UCA'      %Define channel for UCA      
            SimParams.fixedGroupLocs = linspace(-179,180,SimParams.fixedGroups + 1);
            SimParams.fixedGroupLocs = SimParams.fixedGroupLocs(2:2:SimParams.fixedGroups);
            userClustersPerGroup = ceil(length(SimParams.fixedGroupLocs) / SimParams.nGroups);
    end
    
    SimParams.angSpread = inParams.uAngularSpread;
    SimParams.antennasPerGroup = SimParams.nTransmit / SimParams.nGroups; %This is not required
    
    userIndex = 0;
    userIndices = randperm(SimParams.nUsers);
	%%%%%%%% for 45 degree group angle %%%%%%%%%%%%
%     for iGroup = 1:SimParams.nGroups %Defining the group specific terms
%         possibleBeams = SimParams.fixedGroupLocs(((iGroup - 1) * userClustersPerGroup + 5) : (min(iGroup * userClustersPerGroup,length(SimParams.fixedGroupLocs)) - 5)); %Possible number of stat beams to each group
%         SimParams.userLocs{iGroup,1} = possibleBeams(randi([1, length(possibleBeams)],1,SimParams.usersPerGroup(1,iGroup))); %Users location
%         SimParams.groupUserIndices{iGroup,1} = userIndices((userIndex + 1) : sum(SimParams.usersPerGroup(1:iGroup))); %Users that fall into one group
%         SimParams.groupUserBaseAngles{iGroup,1} = possibleBeams; %Start and end angle
%         userIndex = sum(SimParams.usersPerGroup(1:iGroup)); %Users index
%         SimParams.elevAngle{iGroup,1} = linspace(30,60,SimParams.usersPerGroup(1,iGroup)); %Required while considering UCA setup
%     end


	%%%%%%%%%%%for 30 degree group angle %%%%%%%%%%%%
    groupBeamAngles = reshape(SimParams.fixedGroupLocs,[],SimParams.nGroups).';
    for iGroup = 1:SimParams.nGroups %Defining the group specific terms
        possibleBeams = groupBeamAngles(iGroup,42 : end - 42);
        SimParams.userLocs{iGroup,1} = possibleBeams(randi([1, length(possibleBeams)],1,SimParams.usersPerGroup(1,iGroup)));
        SimParams.groupUserIndices{iGroup,1} = userIndices((userIndex + 1) : sum(SimParams.usersPerGroup(1:iGroup))); %Users that fall into one group
        SimParams.groupUserBaseAngles{iGroup,1} = possibleBeams; %Start and end angle
        userIndex = sum(SimParams.usersPerGroup(1:iGroup)); %Users index
        SimParams.elevAngle{iGroup,1} = linspace(30,60,SimParams.usersPerGroup(1,iGroup)); %Required while considering UCA setup
    end
    
    SimParams.groupInfo = struct();
    for iGroup = 1:SimParams.nGroups %Defining struct for one specific group
        SimParams.groupInfo(iGroup).nUsers = SimParams.usersPerGroup(1,iGroup);
        SimParams.groupInfo(iGroup).userLocs = SimParams.userLocs{iGroup,1};
        SimParams.groupInfo(iGroup).elevationAngle = SimParams.elevAngle{iGroup,1}; %UCA specific
        SimParams.groupInfo(iGroup).baseTheta = SimParams.groupUserBaseAngles{iGroup,1};
        SimParams.groupInfo(iGroup).gUserIndices = SimParams.groupUserIndices{iGroup,1};
        SimParams.groupInfo(iGroup).activeBeams = SimParams.gStatBeams;
        SimParams.groupInfo(iGroup).userChannel = zeros(SimParams.nReceive,SimParams.nTransmit,SimParams.groupInfo(iGroup).nUsers); %Initialization
    end
    
    SimParams = outerBeamformerDesign(SimParams); %Find the Outer Precoder B matrix
    SimParams = digitalBeamformerDesign(SimParams); %Find the inner precoder w matrix    
    
    if ~exist(folderName,'dir')
        mkdir(folderName);
    end
    
    fileName = sprintf('%s%stotBeams_%d.mat',folderName,filesep,iBeamPerGroup);
    save(fileName,'SimParams');  
    
    clear SimParams;
    
end

end

