function SimParams = UCA(SimParams)


M = SimParams.nTransmit;
D = sqrt((1-cos(2*pi/M))^2 + sin(2*pi/M)^2);
 SimParams.radius = (SimParams.lambda * 0.5) / D;