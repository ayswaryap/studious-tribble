function B(xlFileName,sheetIDs,varargin)

addpath(genpath(pwd));
warning('off','all')

saveDirectory = sprintf('%s%s%s%s%s%s%s',pwd,filesep,'30degree',filesep,date);

if ~exist(saveDirectory,'dir')
    mkdir(saveDirectory);
end

[~ ,workSheets] = xlsfinfo(xlFileName);
sheetConfig = cell(length(workSheets),1);
copyfile(xlFileName,saveDirectory);

for iSheet = 1:length(sheetIDs)
    
    workBookFolder = sprintf('%s%s%s',saveDirectory,filesep,workSheets{1,sheetIDs(1,iSheet)});
    if ~exist(workBookFolder,'dir')
        mkdir(workBookFolder);
    end
    fileName = sprintf('%s%sGlobalConfig.mat',workBookFolder,filesep);
    
    [xConfig.num, xConfig.txt, xConfig.raw] = xlsread(xlFileName,workSheets{1,sheetIDs(1,iSheet)});
    sheetConfig{iSheet,1} = parseXLFile(xConfig);
    for iColumn = 1:length(sheetConfig{iSheet,1})
        xFields = fieldnames(sheetConfig{iSheet,1}{iColumn,1});
        for iField = 2:length(xFields)
            if ischar(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1})) && contains(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}),'[')
                sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}) = eval(sheetConfig{iSheet,1}{iColumn,1}.(xFields{iField,1}));
            end
        end
    end
    
    if isempty(varargin)
        colIndices = 1:length(sheetConfig{iSheet,1});
    else
        colIndices = varargin{1,1};
    end
   
    save(fileName,'sheetConfig');    
    curWorkSheet = sheetConfig{iSheet,1};
    if isunix
        parfor iColumn = colIndices
            folderName = sprintf('%s%sBookColumn_%d',workBookFolder,filesep,iColumn);
            performTwoLevelPrecoding(curWorkSheet{iColumn,1},folderName);
        end
    else        
        for iColumn = colIndices
            folderName = sprintf('%s%sBookColumn_%d',workBookFolder,filesep,iColumn);
            TLP(curWorkSheet{iColumn,1},folderName);
        end
    end
end

end



