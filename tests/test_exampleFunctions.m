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

                codeOnly = regexprep(lines, "%.*$", "");
                clearMatches = regexp(cellstr(codeOnly), ...
                    '\<clear(?:vars)?\>', 'once');
                testCase.verifyFalse(any(~cellfun('isempty', clearMatches)), ...
                    sprintf('%s must not clear its function workspace.', filePath));
            end
        end

    end
end
