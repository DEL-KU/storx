function varargout = runExample(exampleFile, varargin)
%RUNEXAMPLE Run an example function selected by its file path.
%   RUNEXAMPLE(EXAMPLEFILE) runs an example relative to 00-examples.
%   RUNEXAMPLE(EXAMPLEFILE, NAME, VALUE, ...) forwards name-value inputs to
%   the example. The example folder temporarily takes precedence on the
%   MATLAB path so examples with the same filename remain unambiguous.

if ~(ischar(exampleFile) && isrow(exampleFile)) && ...
        ~(isstring(exampleFile) && isscalar(exampleFile))
    error('storx:examples:InvalidExampleFile', ...
        'exampleFile must be a character vector or string scalar.');
end

examplesRoot = fileparts(mfilename('fullpath'));
exampleFile = char(exampleFile);
if isempty(fileparts(exampleFile)) || ~isfile(exampleFile)
    exampleFile = fullfile(examplesRoot, exampleFile);
end
if isempty(fileparts(exampleFile)) || ~isfile(exampleFile)
    error('storx:examples:ExampleNotFound', ...
        'Example file not found: %s', exampleFile);
end

[ok, attributes] = fileattrib(exampleFile);
if ~ok
    error('storx:examples:ExampleNotFound', ...
        'Example file not found: %s', exampleFile);
end
exampleFile = attributes.Name;
[exampleFolder, functionName, extension] = fileparts(exampleFile);
if ~strcmpi(extension, '.m')
    error('storx:examples:InvalidExampleFile', ...
        'Example file must have a .m extension: %s', exampleFile);
end

originalFolder = pwd;
originalPath = path;
cleanup = onCleanup(@() restoreEnvironment(originalFolder, originalPath));

addpath(exampleFolder, '-begin');
cd(exampleFolder);
rehash path;

resolvedFile = which(functionName);
if isempty(resolvedFile) || ~sameFile(resolvedFile, exampleFile)
    error('storx:examples:ResolutionFailed', ...
        'Could not resolve %s to %s.', functionName, exampleFile);
end

if nargout == 0
    feval(functionName, varargin{:});
else
    [varargout{1:nargout}] = feval(functionName, varargin{:});
end

end

function restoreEnvironment(originalFolder, originalPath)
path(originalPath);
cd(originalFolder);
end

function tf = sameFile(firstFile, secondFile)
if ispc
    tf = strcmpi(firstFile, secondFile);
else
    tf = strcmp(firstFile, secondFile);
end
end
