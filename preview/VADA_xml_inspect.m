%% VADA_XML_INSPECT.m
% Dumps the COMPLETE structure of a .vada.xml and .vada.cfg.xml file
% so you can see exactly what metadata Vevo exports.
%
% Run this on ONE file first to discover which fields contain
% study/series information, then update VADA_organize.m if needed.
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

%% CONFIGURATION
% Point to ONE specific VADA file (without extension)
dataFolder   = 'C:\path\to\VADA_data';  % CHANGE
baseFilename = '';  % CHANGE (no extension, e.g., 'Study_2026-02-04')
vadaExt      = '.vada';

%% INSPECT .vada.xml
xmlPath = fullfile(dataFolder, [baseFilename vadaExt '.xml']);
fprintf('========================================\n');
fprintf('INSPECTING: %s\n', xmlPath);
fprintf('========================================\n\n');

if ~exist(xmlPath, 'file')
    fprintf('FILE NOT FOUND: %s\n', xmlPath);
else
    % Raw text dump of interesting lines
    fprintf('--- RAW TEXT SEARCH (study/series/image/date keywords) ---\n');
    xmlText = fileread(xmlPath);
    lines = strsplit(xmlText, {'\n', '\r'});
    keywords = {'study', 'series', 'image', 'acqui', 'date', 'name', 'label', ...
                'session', 'experiment', 'scan', 'protocol'};
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line), continue; end
        lineLower = lower(line);
        for k = 1:numel(keywords)
            if contains(lineLower, keywords{k})
                fprintf('  Line %4d: %s\n', i, line);
                break;
            end
        end
    end
    
    % Full XML tree dump
    fprintf('\n--- FULL XML TREE ---\n');
    try
        xmlDoc = xmlread(xmlPath);
        dump_xml_node(xmlDoc.getDocumentElement(), 0);
    catch ME
        fprintf('XML parse error: %s\n', ME.message);
    end
end

%% INSPECT .vada.cfg.xml
cfgPath = fullfile(dataFolder, [baseFilename vadaExt '.cfg.xml']);
fprintf('\n\n========================================\n');
fprintf('INSPECTING: %s\n', cfgPath);
fprintf('========================================\n\n');

if ~exist(cfgPath, 'file')
    fprintf('FILE NOT FOUND: %s\n', cfgPath);
else
    fprintf('--- FULL XML TREE ---\n');
    try
        xmlDoc = xmlread(cfgPath);
        dump_xml_node(xmlDoc.getDocumentElement(), 0);
    catch ME
        fprintf('XML parse error: %s\n', ME.message);
    end
end

%% CHECK FOLDER STRUCTURE
fprintf('\n\n========================================\n');
fprintf('FOLDER STRUCTURE around data file\n');
fprintf('========================================\n\n');

% Show what's in the data folder
fprintf('Contents of %s:\n', dataFolder);
d = dir(dataFolder);
for i = 1:numel(d)
    if d(i).isdir
        fprintf('  [DIR]  %s\n', d(i).name);
    else
        fprintf('  %6.1f MB  %s\n', d(i).bytes/1e6, d(i).name);
    end
end

% Check parent folder (Vevo sometimes organizes study > series > files)
parentDir = fileparts(dataFolder);
if ~isempty(parentDir)
    fprintf('\nParent (%s):\n', parentDir);
    d = dir(parentDir);
    for i = 1:min(20, numel(d))  % Cap at 20 entries
        if d(i).isdir
            fprintf('  [DIR]  %s\n', d(i).name);
        else
            fprintf('  %6.1f MB  %s\n', d(i).bytes/1e6, d(i).name);
        end
    end
end

fprintf('\n=== Inspection complete ===\n');
fprintf('Look for study/series fields in the XML tree dump above.\n');
fprintf('If found, VADA_organize.m should pick them up automatically.\n');
fprintf('If NOT found, the data may use folder structure or VevoLab DB only.\n');

%% ========================================================================
function dump_xml_node(node, depth)
% Recursively dump an XML node tree with indentation

indent = repmat('  ', 1, depth);
nodeName = char(node.getNodeName);

% Print node with attributes
attrStr = '';
if node.hasAttributes
    attrs = node.getAttributes;
    for i = 0:attrs.getLength-1
        attr = attrs.item(i);
        attrStr = [attrStr sprintf(' %s="%s"', char(attr.getName), char(attr.getValue))]; %#ok
    end
end

% Get text content (if leaf node)
textContent = '';
if node.hasChildNodes
    children = node.getChildNodes;
    if children.getLength == 1 && children.item(0).getNodeType == 3  % TEXT_NODE
        textContent = strtrim(char(children.item(0).getTextContent));
    end
end

if ~isempty(textContent)
    fprintf('%s<%s%s> %s\n', indent, nodeName, attrStr, textContent);
else
    fprintf('%s<%s%s>\n', indent, nodeName, attrStr);
end

% Recurse into children (skip text nodes)
if node.hasChildNodes
    children = node.getChildNodes;
    for i = 0:children.getLength-1
        child = children.item(i);
        if child.getNodeType == 1  % ELEMENT_NODE
            dump_xml_node(child, depth + 1);
        end
    end
end
end
