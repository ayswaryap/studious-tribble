function plotResults(matFile,plotType)

if ~iscell(matFile)
    matFile = {matFile};
end

outLegend = {};
xFiles = cell(length(matFile),1);
markerTypes = {'*','o','d','s','v','h','^'};

for iFile = 1:length(matFile)
    
    load(matFile{iFile});
   
    iResponse = 0;
    for iColumn = 1:length(xParams)
        for iRow = 1:length(xParams{iColumn,1})
            iResponse = iResponse + 1;
            xFiles{iFile,1}(iResponse) = xParams{iColumn,1}{iRow,1};
        end        
    end   
   
    switch plotType
        case 'Beams'
            outLegend = plotThrptVsBeams(xFiles{iFile,1}, outLegend, markerTypes{1,iFile});
    end
    
    legend(outLegend);

end

end

function [outLegend] = plotThrptVsBeams(SimParams, outLegend, markerType)
    
    sumRate = zeros(1,length(SimParams));
    beamCount = zeros(1,length(SimParams));
    for iConfig = 1:length(SimParams)
        xRate = SimParams(iConfig).groupSumRate.srate(1:SimParams(iConfig).groupSumRate.tsca);
        sumRate(1,iConfig) = mean(max(xRate,[],2),1);
        beamCount(1,iConfig) = SimParams(iConfig).gStatBeams * SimParams(iConfig).nGroups;
    end
    
    plot(beamCount,sumRate,'Marker',markerType);
    outLegend{1,length(outLegend) + 1} = sprintf('%s-%s with groups and beams <%d, %d>',SimParams(1).statBeamType,SimParams(1).innerPrecoder,SimParams(1).nGroups,SimParams(1).gStatBeams * SimParams(1).nGroups);
    hold all;
    
end


