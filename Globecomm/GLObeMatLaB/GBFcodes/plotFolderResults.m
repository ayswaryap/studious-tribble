
function legendName = plotFolderResults(varargin)

switch size(varargin,2)
    case 0
        parentFolder = 'Results\14-Feb-2018'; %'AnttiProjects\GroupBeamforming\TJ6OSC~1\05-Jan-2018';
        workSheetNos = [];
    case 1
        parentFolder = varargin{1,1};
        workSheetNos = [];
    case 2
        parentFolder = varargin{1,1};
        workSheetNos = varargin{1,2};
end

lsDir = dir(parentFolder);
xlFileID = arrayfun(@(x)(~isempty(strfind(lsDir(x).name,'.xlsx')) && isempty(strfind(lsDir(x).name,'~$'))),1:length(lsDir));

if any(xlFileID)
    xlFileName = sprintf('%s%s%s',lsDir(xlFileID).folder,filesep,lsDir(xlFileID).name);
    [~ ,workSheets] = xlsfinfo(xlFileName);
    if ~isempty(workSheetNos)
        workSheets = workSheets(workSheetNos);
    end
%     [xConfig.num, xConfig.txt, xConfig.raw] = xlsread(xlFileName,workSheets{1,sheetIDs(1,iSheet)});
%     sheetConfig{iSheet,1} = parseXLFile(xConfig);
end

legendName = {};
for iWorkBook = 1:length(lsDir)
    if any(strcmpi(lsDir(iWorkBook).name,workSheets))
        lsCols = dir(sprintf('%s%s%s',lsDir(iWorkBook).folder,filesep,lsDir(iWorkBook).name));
        for iCol = 1:length(lsCols)
            if ~any(strcmpi(lsCols(iCol).name,{'.','..'})) && lsCols(iCol).isdir
                lsBeams = dir(sprintf('%s%s%s',lsCols(iCol).folder,filesep,lsCols(iCol).name));
                beamFiles = lsBeams(arrayfun(@(x)(~x.isdir),lsBeams));
                legendName = displayPlots(beamFiles,lsDir(iWorkBook).name,legendName);
            end
        end
    end
end

display(legendName);

end

function legendName = displayPlots(beamFiles,workSheetName,legendName)

if isempty(beamFiles)
    return;
end

persistent markerIndex;

if isempty(markerIndex)
    markerIndex = 1;
else
    markerIndex = markerIndex + 1;
end

markerType = {'+','*','o','.','p','v','d','s','>','^','x'};

legendIndex = length(legendName) + 1;
fprintf('Processing worksheet : %s \n',workSheetName);

beamPlotX = zeros(1,length(beamFiles));
beamPlotY = zeros(1,length(beamFiles));
beamPlotZ = zeros(1,length(beamFiles));

for iBeam = 1:length(beamFiles)
    
    clear SimParams;
    load(sprintf('%s%s%s',beamFiles(iBeam).folder,filesep,beamFiles(iBeam).name));
    if iBeam == 1
        legendName{1,legendIndex} = sprintf('OBF - %s IBF - %s : [ G / U_g ] : [%d / %d] - %d dB',SimParams.statBeamType,SimParams.innerPrecoder,SimParams.nGroups,SimParams.groupInfo(1).nUsers,db(SimParams.txPower,'power'));
    end
    
    validFields = logical(arrayfun(@(x)(SimParams.groupSumRate.isSucceded(x)),1:SimParams.nMontRuns));
    if any(validFields)
        userSumRate = max(SimParams.groupSumRate.srate(validFields,:),[],2);
        beamPlotY(1,iBeam) = mean(userSumRate);
    else
        beamPlotY(1,iBeam) = nan;
    end
        
    fprintf('Number of valid samples : %d for :: %s \n',sum(validFields),legendName{1,legendIndex});
    beamPlotX(1,iBeam) = SimParams.totStatBeams;%SimParams.beamsPerGroup;
    if isfield(SimParams.groupSumRate,'stime')
        beamPlotZ(1,iBeam) = mean(SimParams.groupSumRate.stime(validFields,1),1);
    end
    
end

hold all;
figure(1);
if size(beamPlotX) == 1
    plot((1:beamPlotX(~isnan(beamPlotY))),repmat(beamPlotY(~isnan(beamPlotY)),1,beamPlotX(~isnan(beamPlotY))),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
else
    plot(beamPlotX(~isnan(beamPlotY)),beamPlotY(~isnan(beamPlotY)),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
end

% figure(2);
% plot(beamPlotX(~isnan(beamPlotY)),beamPlotZ(~isnan(beamPlotZ)),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
% 
fprintf('Displaying performance of [%s] \n',legendName{1,legendIndex});

end




















%%%%%%%%%%%%%%%Use to plot the complexity


% 
% function legendName = plotFolderResults(varargin)
% 
% switch size(varargin,2)
%     case 0
%         parentFolder = 'Results\14-Feb-2018'; %'AnttiProjects\GroupBeamforming\TJ6OSC~1\05-Jan-2018';
%         workSheetNos = [];
%     case 1
%         parentFolder = varargin{1,1};
%         workSheetNos = [];
%     case 2
%         parentFolder = varargin{1,1};
%         workSheetNos = varargin{1,2};
% end
% 
% lsDir = dir(parentFolder);
% xlFileID = arrayfun(@(x)(~isempty(strfind(lsDir(x).name,'.xlsx')) && isempty(strfind(lsDir(x).name,'~$'))),1:length(lsDir));
% 
% if any(xlFileID)
%     xlFileName = sprintf('%s%s%s',lsDir(xlFileID).folder,filesep,lsDir(xlFileID).name);
%     [~ ,workSheets] = xlsfinfo(xlFileName);
%     if ~isempty(workSheetNos)
%         workSheets = workSheets(workSheetNos);
%     end
%     [xConfig.num, xConfig.txt, xConfig.raw] = xlsread(xlFileName,workSheets{1,sheetIDs(1,iSheet)});
%     sheetConfig{iSheet,1} = parseXLFile(xConfig);
% end
% 
% legendName = {};
% for iWorkBook = 1:length(lsDir)
%     if any(strcmpi(lsDir(iWorkBook).name,workSheets))
%         lsCols = dir(sprintf('%s%s%s',lsDir(iWorkBook).folder,filesep,lsDir(iWorkBook).name));
%         for iCol = 1:length(lsCols)
%             if ~any(strcmpi(lsCols(iCol).name,{'.','..'})) && lsCols(iCol).isdir
%                 lsBeams = dir(sprintf('%s%s%s',lsCols(iCol).folder,filesep,lsCols(iCol).name));
%                 beamFiles = lsBeams(arrayfun(@(x)(~x.isdir),lsBeams));
%                 legendName = displayPlots(beamFiles,lsDir(iWorkBook).name,legendName);
%             end
%         end
%     end
% end
% 
% display(legendName);
% 
% end
% 
% function legendName = displayPlots(beamFiles,workSheetName,legendName)
% 
% if isempty(beamFiles)
%     return;
% end
% 
% persistent markerIndex;
% 
% if isempty(markerIndex)
%     markerIndex = 1;
% else
%     markerIndex = markerIndex + 1;
% end
% 
% markerType = {'+','*','o','.','p','v','d','s','>','^','x'};
% 
% legendIndex = length(legendName) + 1;
% fprintf('Processing worksheet : %s \n',workSheetName);
% 
% beamPlotX = zeros(1,length(beamFiles));
% beamPlotY = zeros(1,length(beamFiles));
% beamPlotZ = zeros(1,length(beamFiles));
% 
% for iBeam = 1:length(beamFiles)
%     
%     clear SimParams;
%     load(sprintf('%s%s%s',beamFiles(iBeam).folder,filesep,beamFiles(iBeam).name));
%     if iBeam == 1
%         legendName{1,legendIndex} = sprintf('%s-%d-%s-%s',workSheetName,SimParams.usersPerGroup(1),SimParams.statBeamType,SimParams.innerPrecoder);
%     end
%     
%     validFields = logical(arrayfun(@(x)(SimParams.groupSumRate.isSucceded(x)),1:SimParams.nMontRuns));
%     if any(validFields)
%         userSumRate = max(SimParams.groupSumRate.srate(validFields,:),[],2);
%         beamPlotY(1,iBeam) = mean(userSumRate);
%     else
%         beamPlotY(1,iBeam) = nan;
%     end
%         
%     fprintf('Number of valid samples : %d for :: %s \n',sum(validFields),legendName{1,legendIndex});
%     beamPlotX(1,iBeam) = SimParams.beamsPerGroup;
%     if isfield(SimParams.groupSumRate,'stime')
%         beamPlotZ(1,iBeam) = mean(SimParams.groupSumRate.stime(validFields,1),1);
%     end
%     
% end
% 
% hold all;
% figure(1);
% if size(beamPlotX) == 1
%     plot((1:beamPlotX(~isnan(beamPlotY))),repmat(beamPlotY(~isnan(beamPlotY)),1,beamPlotX(~isnan(beamPlotY))),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
% else
%     plot(beamPlotX(~isnan(beamPlotY)),beamPlotY(~isnan(beamPlotY)),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
% end
% 
% figure(2);
% plot(beamPlotX(~isnan(beamPlotY)),beamPlotZ(~isnan(beamPlotZ)),'Marker',markerType{1,mod(markerIndex,length(markerType)) + 1});
% 
% fprintf('Displaying performance of [%s] \n',legendName{1,legendIndex});
% 
% end
% 
% 
% 
% 
% 
