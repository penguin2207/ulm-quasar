%% VADA_ORGANIZE.m
% Organize VADA export files into folders by Study and Series.
%
% Reads from .vada.xml parameter list (format confirmed via VADA_xml_inspect):
%   <parameter name="Study-Name" value="Study (2026-02-04 16:56:50)"/>
%   <parameter name="Series-Name" value="Standard"/>
%
% Output structure:
%   outputRoot/
%     Study_(2026-02-04_16-56-50)/
%       Standard/
%         block1.vada
%         block1.vada.xml
%         block1.vada.cfg.xml
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

%% CONFIGURATION
inputFolder = 'C:\path\to\VADA_data';  % CHANGE
outputRoot  = fullfile(inputFolder, 'organized');
dryRun      = true;    % true = preview only, false = execute
useMove     = false;    % false = copy (safer), true = move
vadaExt     = '.vada';

%% DISCOVER FILES
fprintf('=== VADA File Organizer ===\n');
fprintf('Input:  %s\n', inputFolder);
fprintf('Output: %s\n', outputRoot);
fprintf('Mode:   %s\n\n', iff(dryRun, 'DRY RUN (preview)', iff(useMove, 'MOVE', 'COPY')));

% Find .vada files (top level + one level of subdirs, skip output dirs)
skipDirs = {'organized', 'LAT_ULM_Results'};
vadaFiles = dir(fullfile(inputFolder, ['*' vadaExt]));
subDirs = dir(inputFolder);
subDirs = subDirs([subDirs.isdir] & ~startsWith({subDirs.name}, '.') & ...
                  ~ismember({subDirs.name}, skipDirs));
for i = 1:numel(subDirs)
    vadaFiles = [vadaFiles; dir(fullfile(inputFolder, subDirs(i).name, ['*' vadaExt]))]; %#ok
end

if isempty(vadaFiles)
    fprintf('No %s files found.\n', vadaExt); return;
end
fprintf('Found %d .vada file(s)\n\n', numel(vadaFiles));

%% EXTRACT METADATA
fprintf('--- Reading metadata ---\n\n');

nFiles = numel(vadaFiles);
study    = cell(nFiles, 1);
series   = cell(nFiles, 1);
seqName  = cell(nFiles, 1);
sos      = cell(nFiles, 1);
voltage  = cell(nFiles, 1);
folders  = cell(nFiles, 1);
baseNames = cell(nFiles, 1);
sizes    = zeros(nFiles, 1);

for i = 1:nFiles
    folders{i}   = vadaFiles(i).folder;
    baseNames{i} = vadaFiles(i).name(1:end-numel(vadaExt));
    sizes(i)     = vadaFiles(i).bytes / 1e9;
    
    % Read parameters from .vada.xml
    xmlPath = fullfile(folders{i}, [baseNames{i} vadaExt '.xml']);
    params = read_vada_xml_params(xmlPath);
    
    study{i}   = get_param(params, 'Study-Name', '');
    series{i}  = get_param(params, 'Series-Name', '');
    seqName{i} = get_param(params, 'Vada-Mode/User-Pulse-Sequence-Name', '');
    sos{i}     = get_param(params, 'Vada-Mode/Speed-Of-Sound-Media', '');
    
    vLo = get_param(params, 'Vada-Mode/Voltage-Rail-Low', '');
    vHi = get_param(params, 'Vada-Mode/Voltage-Rail-High', '');
    if ~isempty(vLo)
        voltage{i} = sprintf('%s-%s%%', vLo, vHi);
    end
    
    fprintf('  [%2d] %s  (%.1f GB)\n', i, vadaFiles(i).name, sizes(i));
    fprintf('       Study:  %s\n', study{i});
    fprintf('       Series: %s | Seq: %s | SoS: %s | V: %s\n', ...
        series{i}, seqName{i}, sos{i}, voltage{i});
end

%% GROUP BY STUDY + SERIES
fprintf('\n--- Grouping ---\n\n');

groupKeys = cell(nFiles, 1);
for i = 1:nFiles
    s = study{i};  if isempty(s), s = 'UnknownStudy'; end
    r = series{i}; if isempty(r), r = 'UnknownSeries'; end
    groupKeys{i} = [s '||' r];
end

uniqueGroups = unique(groupKeys, 'stable');
fprintf('Found %d group(s):\n\n', numel(uniqueGroups));

for g = 1:numel(uniqueGroups)
    members = find(strcmp(groupKeys, uniqueGroups{g}));
    parts = strsplit(uniqueGroups{g}, '||');
    studyName = parts{1}; seriesName = parts{2};
    totalGB = sum(sizes(members));
    
    fprintf('  Group %d: "%s" / "%s"\n', g, studyName, seriesName);
    fprintf('    %d block(s), %.1f GB total\n', numel(members), totalGB);
    for m = members'
        fprintf('      %s\n', vadaFiles(m).name);
    end
    fprintf('\n');
end

%% ORGANIZE
fprintf('--- Organization plan ---\n\n');

nCopied = 0;
for g = 1:numel(uniqueGroups)
    members = find(strcmp(groupKeys, uniqueGroups{g}));
    parts = strsplit(uniqueGroups{g}, '||');
    
    studyFolder  = safe_name(parts{1});
    seriesFolder = safe_name(parts{2});
    destDir = fullfile(outputRoot, studyFolder, seriesFolder);
    
    fprintf('  %s/%s/\n', studyFolder, seriesFolder);
    
    for m = members'
        companions = find_companions(folders{m}, baseNames{m}, vadaExt);
        
        for c = 1:numel(companions)
            srcPath = fullfile(folders{m}, companions{c});
            dstPath = fullfile(destDir, companions{c});
            
            fprintf('    %s\n', companions{c});
            
            if ~dryRun
                if ~exist(destDir, 'dir'), mkdir(destDir); end
                if useMove
                    movefile(srcPath, dstPath);
                else
                    copyfile(srcPath, dstPath);
                end
                nCopied = nCopied + 1;
            end
        end
    end
    fprintf('\n');
end

if dryRun
    fprintf('*** DRY RUN complete. Set dryRun=false to execute. ***\n');
else
    fprintf('*** Done: %d files %s to %s ***\n', ...
        nCopied, iff(useMove, 'moved', 'copied'), outputRoot);
end

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function params = read_vada_xml_params(xmlPath)
% Read all <parameter name="..." value="..."/> entries from a .vada.xml file.
% Returns a struct mapping name -> value (as strings).
params = struct();

if ~exist(xmlPath, 'file'), return; end

text = fileread(xmlPath);

% Match: <parameter name="X" value="Y" .../>
tokens = regexp(text, '<parameter\s+name="([^"]+)"\s+value="([^"]*)"', 'tokens');

for i = 1:numel(tokens)
    params.(tokens{i}{1}) = tokens{i}{2};
end
end


function val = get_param(params, name, default)
% Get a parameter value, handling the hyphenated names from Vevo XML.
if isstruct(params)
    % Try exact field name first (with hyphens replaced by underscores for struct)
    fn = strrep(strrep(name, '-', '_'), '/', '_');
    if isfield(params, fn)
        val = params.(fn);
        return;
    end
    % Try original name as a dynamic field match
    fnames = fieldnames(params);
    for i = 1:numel(fnames)
        if strcmp(strrep(fnames{i}, '_', '-'), name) || ...
           strcmpi(fnames{i}, strrep(strrep(name, '-', ''), '/', ''))
            val = params.(fnames{i});
            return;
        end
    end
end
val = default;
end


function companions = find_companions(folder, baseName, vadaExt)
% Find all files belonging to a VADA block.
companions = {};
exts = {vadaExt, [vadaExt '.xml'], [vadaExt '.cfg.xml']};
for e = 1:numel(exts)
    f = [baseName exts{e}];
    if exist(fullfile(folder, f), 'file')
        companions{end+1} = f; %#ok
    end
end
% Any other companion files with same base name
others = dir(fullfile(folder, [baseName '*']));
for f = 1:numel(others)
    if ~ismember(others(f).name, companions) && ~others(f).isdir
        companions{end+1} = others(f).name; %#ok
    end
end
end


function s = safe_name(name)
% Make a filesystem-safe folder name.
s = strrep(name, ':', '-');
s = strrep(s, '/', '_');
s = strrep(s, '\', '_');
s = strrep(s, '"', '');
s = strrep(s, '<', '');
s = strrep(s, '>', '');
s = strrep(s, '|', '_');
s = strrep(s, '?', '');
s = strrep(s, '*', '');
s = strtrim(s);
if isempty(s), s = 'unnamed'; end
end


function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end
