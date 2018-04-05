function twoLevelBeamformerDesign(xlFileName,sheetIDs)


addpath(genpath(pwd));
warning('off','all')

[~ ,workSheets] = xlsfinfo(xlFileName);
sheetConfig = cell(length(workSheets),1);

for iSheet = 1:length(sheetIDs)
    fileName = sprintf('Output_%s.mat',workSheets{1,sheetIDs(1,iSheet)});
    [xConfig.num, xConfig.txt, xConfig.raw] = xlsread(xlFileName,workSheets{1,sheetIDs(1,iSheet)});
    sheetConfig{iSheet,1} = parseXLFile(xConfig);
    for iColumn = 1:length(sheetConfig{iSheet,1})
        xFields = fieldnames(sheetConfig{iSheet,1}{iColumn,1});
        for iField = 2:length(xFields)
            if ischar(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1})) && ~isempty(strfind(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}),'['))
                sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}) = eval(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}));
            end
        end
    end
    
    xParams = cell(length(sheetConfig{iSheet,1}),1);
    for iColumn = 1:length(sheetConfig{iSheet,1}) % parfor
        xParams{iColumn,1} = performTwoLevelPrecoding(sheetConfig{iSheet,1}{iColumn,1});
        save(fileName,'xParams','sheetConfig');
    end
end

save(fileName,'xParams','sheetConfig');

end

