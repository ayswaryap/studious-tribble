function plotBeamPerformance(SimParams,sBeamM)

nResols = 1024;
xTheta = linspace(-179,180,nResols);
plotGains = zeros(nResols,2,size(sBeamM,2));

for iTheta = 1:nResols
    for iBeam = 1:size(sBeamM,2)
        cWBeamformer = exp(-sqrt(-1) * pi * sind(xTheta(1,iTheta)) * linspace(0,SimParams.nTransmit - 1,SimParams.nTransmit)) / sqrt(SimParams.nTransmit);
        plotGains(iTheta,2,iBeam) = cWBeamformer * sBeamM(:,iBeam);
        plotGains(iTheta,1,iBeam) = xTheta(1,iTheta);
    end
end

for iBeam = 1:size(sBeamM,2)
    plot(plotGains(:,1,iBeam),abs(plotGains(:,2,iBeam)),':');hold all;
end

end

