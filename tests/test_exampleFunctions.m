classdef test_exampleFunctions < matlab.unittest.TestCase
    %% Source-level checks for the chapter examples

    methods (Test)

        function examplesUseFunctionWorkspaces(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            repositoryRoot = fileparts(testFolder);
            examplesRoot = fullfile(repositoryRoot, '00-examples');
            files = dir(fullfile(examplesRoot, '**', '*.m'));

            % Root-level files are utilities; nested files are examples.
            files = files(~strcmp({files.folder}, examplesRoot));
            testCase.assertNotEmpty(files, 'No example files were found.');

            for fileIndex = 1:numel(files)
                file = files(fileIndex);
                filePath = fullfile(file.folder, file.name);
                source = string(fileread(filePath));
                lines = splitlines(source);
                trimmedLines = strtrim(lines);

                functionIndex = find(startsWith(trimmedLines, "function "), 1);
                testCase.verifyNotEmpty(functionIndex, sprintf( ...
                    '%s must declare a primary function.', filePath));
                if isempty(functionIndex)
                    continue
                end

                [~, expectedName] = fileparts(file.name);
                testCase.verifyTrue(contains(trimmedLines(functionIndex), ...
                    expectedName + "("), sprintf( ...
                    '%s primary function must match its filename.', filePath));

                remainingIndices = (functionIndex + 1):numel(trimmedLines);
                isCode = strlength(trimmedLines(remainingIndices)) > 0 & ...
                    ~startsWith(trimmedLines(remainingIndices), "%");
                firstCodeIndex = remainingIndices(find(isCode, 1));
                testCase.verifyTrue(startsWith( ...
                    trimmedLines(firstCodeIndex), "arguments"), sprintf( ...
                    '%s must put an arguments block before executable code.', ...
                    filePath));

                argumentsEndOffset = find( ...
                    trimmedLines((firstCodeIndex + 1):end) == "end", 1);
                testCase.verifyNotEmpty(argumentsEndOffset, sprintf( ...
                    '%s arguments block must have a matching end.', filePath));
                if isempty(argumentsEndOffset)
                    continue
                end

                argumentsEndIndex = firstCodeIndex + argumentsEndOffset;
                remainingIndices = (argumentsEndIndex + 1):numel(trimmedLines);
                isCode = strlength(trimmedLines(remainingIndices)) > 0 & ...
                    ~startsWith(trimmedLines(remainingIndices), "%");
                configureIndex = remainingIndices(find(isCode, 1));
                testCase.verifyEqual(trimmedLines(configureIndex), ...
                    "configureGraphics();", sprintf( ...
                    ['%s must call configureGraphics() immediately after its ' ...
                    'arguments block.'], filePath));

                codeOnly = regexprep(lines, "%.*$", "");
                clearMatches = regexp(cellstr(codeOnly), ...
                    '\<clear(?:vars)?\>', 'once');
                testCase.verifyFalse(any(~cellfun('isempty', clearMatches)), ...
                    sprintf('%s must not clear its function workspace.', filePath));

                clcMatches = regexp(cellstr(codeOnly), '\<clc\>', 'once');
                testCase.verifyFalse(any(~cellfun('isempty', clcMatches)), ...
                    sprintf('%s must not clear the command window.', filePath));

                configureMatches = regexp(cellstr(codeOnly), ...
                    '^\s*configureGraphics\(\);\s*$', 'once');
                testCase.verifyEqual(sum(~cellfun('isempty', configureMatches)), ...
                    1, sprintf('%s must configure graphics exactly once.', filePath));
            end
        end

        function combinedFiguresAreOptional(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            repositoryRoot = fileparts(testFolder);
            examplesRoot = fullfile(repositoryRoot, '00-examples');
            files = dir(fullfile(examplesRoot, '**', '*.m'));

            % Root-level files are utilities; nested files are examples.
            files = files(~strcmp({files.folder}, examplesRoot));
            numExportExamples = 0;

            for fileIndex = 1:numel(files)
                file = files(fileIndex);
                filePath = fullfile(file.folder, file.name);
                lines = splitlines(string(fileread(filePath)));
                codeOnly = regexprep(lines, "%.*$", "");
                trimmedCode = strtrim(codeOnly);

                exportDeclarations = regexp(cellstr(codeOnly), ...
                    '^\s*options\.exportImages\>', 'once');
                if ~any(~cellfun('isempty', exportDeclarations))
                    continue
                end
                numExportExamples = numExportExamples + 1;

                combineDeclarations = find(trimmedCode == ...
                    "options.combineFigure (1,1) logical = false");
                testCase.verifyNumElements(combineDeclarations, 1, sprintf( ...
                    ['%s must declare exactly one logical combineFigure ' ...
                    'option that defaults to false.'], filePath));

                combineAssignments = find(trimmedCode == ...
                    "combineFigure = options.combineFigure;");
                testCase.verifyNumElements(combineAssignments, 1, sprintf( ...
                    ['%s must assign options.combineFigure to the local ' ...
                    'combineFigure variable exactly once.'], filePath));

                combineMatches = regexp(cellstr(codeOnly), ...
                    '\<combineFigures\s*\(', 'once');
                combineIndices = find(~cellfun('isempty', combineMatches));
                testCase.verifyNumElements(combineIndices, 1, sprintf( ...
                    '%s must call combineFigures exactly once.', filePath));

                guardIndices = find(trimmedCode == "if combineFigure");
                testCase.verifyNumElements(guardIndices, 1, sprintf( ...
                    ['%s must guard its combined title, figure call, and save ' ...
                    'with exactly one if combineFigure block.'], filePath));
                if numel(guardIndices) ~= 1
                    continue
                end

                guardIndex = guardIndices(1);
                guardEndIndex = test_exampleFunctions.findBlockEnd( ...
                    trimmedCode, guardIndex);
                testCase.verifyNotEmpty(guardEndIndex, sprintf( ...
                    '%s combineFigure guard must have a matching end.', filePath));
                if isempty(guardEndIndex)
                    continue
                end

                guardedIndices = (guardIndex + 1):(guardEndIndex - 1);
                titleMatches = regexp(cellstr(codeOnly), ...
                    '^\s*ex_title\s*=', 'once');
                titleIndices = find(~cellfun('isempty', titleMatches));
                testCase.verifyNumElements(titleIndices, 1, sprintf( ...
                    '%s must build ex_title exactly once.', filePath));
                guardedTitleIndices = intersect(titleIndices, guardedIndices);
                testCase.verifyNumElements(guardedTitleIndices, 1, sprintf( ...
                    '%s must build ex_title inside its combineFigure guard.', ...
                    filePath));

                guardedCombineIndices = intersect( ...
                    combineIndices, guardedIndices);
                testCase.verifyNumElements(guardedCombineIndices, 1, sprintf( ...
                    '%s must call combineFigures inside its combineFigure guard.', ...
                    filePath));

                nestedExportIndices = guardedIndices( ...
                    trimmedCode(guardedIndices) == "if exportImages");
                testCase.verifyNumElements(nestedExportIndices, 1, sprintf( ...
                    ['%s must guard the post-combine save with if exportImages ' ...
                    'inside its combineFigure guard.'], filePath));
                if numel(nestedExportIndices) ~= 1
                    continue
                end

                exportGuardIndex = nestedExportIndices(1);
                exportGuardEndIndex = test_exampleFunctions.findBlockEnd( ...
                    trimmedCode, exportGuardIndex);
                testCase.verifyNotEmpty(exportGuardEndIndex, sprintf( ...
                    ['%s post-combine exportImages guard must have a matching ' ...
                    'end.'], filePath));
                if isempty(exportGuardEndIndex)
                    continue
                end
                testCase.verifyLessThan(exportGuardEndIndex, guardEndIndex, ...
                    sprintf(['%s post-combine exportImages guard must close ' ...
                    'inside its combineFigure guard.'], filePath));

                saveMatches = regexp(cellstr(codeOnly), ...
                    '^\s*saveAll\s*\(\s*folder\s*\)', 'once');
                saveIndices = find(~cellfun('isempty', saveMatches));
                savesInCombineGuard = intersect(saveIndices, guardedIndices);
                testCase.verifyNumElements(savesInCombineGuard, 1, sprintf( ...
                    ['%s combineFigure guard must contain exactly one ' ...
                    'saveAll(folder) call.'], filePath));
                exportedIndices = (exportGuardIndex + 1): ...
                    (exportGuardEndIndex - 1);
                guardedSaveIndices = intersect(saveIndices, exportedIndices);
                testCase.verifyNumElements(guardedSaveIndices, 1, sprintf( ...
                    ['%s must save the combined figure exactly once inside ' ...
                    'the nested exportImages guard.'], filePath));

                if isscalar(guardedTitleIndices) && ...
                        isscalar(guardedCombineIndices) && ...
                        isscalar(guardedSaveIndices)
                    testCase.verifyLessThan(guardedTitleIndices, ...
                        guardedCombineIndices, sprintf( ...
                        '%s must build ex_title before combining figures.', ...
                        filePath));
                    testCase.verifyLessThan(guardedCombineIndices, ...
                        exportGuardIndex, sprintf( ...
                        ['%s must combine figures before testing whether to ' ...
                        'save them.'], filePath));
                    testCase.verifyLessThan(guardedCombineIndices, ...
                        guardedSaveIndices, sprintf( ...
                        '%s must save after combining figures.', filePath));
                end

                if isscalar(combineAssignments)
                    testCase.verifyLessThan(combineAssignments, guardIndex, ...
                        sprintf(['%s must assign combineFigure before testing ' ...
                        'its guard.'], filePath));
                end
            end

            testCase.verifyEqual(numExportExamples, 83, ...
                'Expected 83 official examples that declare exportImages.');
        end

        function exportFoldersCoverSupportedFormats(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            repositoryRoot = fileparts(testFolder);
            officialRoot = fullfile(repositoryRoot, '00-examples');
            legacyRoots = {
                fullfile(repositoryRoot, 'extras', 'examples_paramOpt')
                fullfile(repositoryRoot, 'extras', 'triFEA', 'examples_fea')
                };

            officialFiles = dir(fullfile(officialRoot, '**', '*.m'));
            officialFiles = officialFiles(~strcmp( ...
                {officialFiles.folder}, officialRoot));
            filePaths = string(fullfile( ...
                {officialFiles.folder}, {officialFiles.name}))';

            for rootIndex = 1:numel(legacyRoots)
                legacyFiles = dir(fullfile(legacyRoots{rootIndex}, '**', '*.m'));
                filePaths = [filePaths; string(fullfile( ... %#ok<AGROW>
                    {legacyFiles.folder}, {legacyFiles.name}))'];
            end

            exportNames = ["exportImages", "exportGIF", "exportSTL"];
            for fileIndex = 1:numel(filePaths)
                filePath = char(filePaths(fileIndex));
                lines = splitlines(string(fileread(filePath)));
                codeOnly = regexprep(lines, "%.*$", "");

                declaredExports = false(size(exportNames));
                for exportIndex = 1:numel(exportNames)
                    declarationPattern = ['^\s*(?:options\.)?' ...
                        char(exportNames(exportIndex)) '\>[^=]*='];
                    declarations = regexp(cellstr(codeOnly), ...
                        declarationPattern, 'once');
                    declaredExports(exportIndex) = any( ...
                        ~cellfun('isempty', declarations));
                end

                mkdirMatches = regexp(cellstr(codeOnly), ...
                    '^\s*mkdir\s*\(\s*folder\s*\)\s*;?\s*$', 'once');
                mkdirIndices = find(~cellfun('isempty', mkdirMatches));

                if ~any(declaredExports)
                    testCase.verifyEmpty(mkdirIndices, sprintf( ...
                        ['%s creates an export folder without declaring an ' ...
                        'export option.'], filePath));
                    continue
                end

                testCase.verifyNumElements(mkdirIndices, 1, sprintf( ...
                    '%s must create exactly one export folder.', filePath));
                if numel(mkdirIndices) ~= 1
                    continue
                end

                mkdirIndex = mkdirIndices(1);
                precedingLines = strtrim(codeOnly(1:mkdirIndex - 1));
                guardIndex = find(startsWith(precedingLines, "if export"), ...
                    1, 'last');
                testCase.verifyNotEmpty(guardIndex, sprintf( ...
                    '%s must guard mkdir(folder) with its export options.', ...
                    filePath));
                if isempty(guardIndex)
                    continue
                end

                expectedGuard = "if " + strjoin( ...
                    exportNames(declaredExports), " || ");
                testCase.verifyEqual(precedingLines(guardIndex), ...
                    expectedGuard, sprintf( ...
                    ['%s export-folder guard must include every supported ' ...
                    'format in canonical order.'], filePath));
            end
        end

    end

    methods (Static, Access = private)

        function blockEndIndex = findBlockEnd(lines, blockStartIndex)
            blockStartPattern = ...
                '^(?:if|for|parfor|while|switch|try|spmd)\>';
            blockEndIndex = [];
            blockDepth = 1;

            for lineIndex = (blockStartIndex + 1):numel(lines)
                line = char(lines(lineIndex));
                if ~isempty(regexp(line, blockStartPattern, 'once'))
                    blockDepth = blockDepth + 1;
                elseif ~isempty(regexp(line, '^end\s*;?$', 'once'))
                    blockDepth = blockDepth - 1;
                    if blockDepth == 0
                        blockEndIndex = lineIndex;
                        return
                    end
                end
            end
        end

    end
end
