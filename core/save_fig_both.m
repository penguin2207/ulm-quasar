function save_fig_both(fig, pngPath, varargin)
%SAVE_FIG_BOTH  Export figure as both PNG (exportgraphics) and .fig (savefig).
%
%   save_fig_both(fig, pngPath, 'Resolution', 150)
%
%  Derives the .fig path by replacing the trailing '.png' with '.fig'.
%  exportgraphics extra args (e.g. 'Resolution', 150) are forwarded.
%  .fig save failures are warnings, not errors, so the batch continues.

figPath = regexprep(pngPath, '\.png$', '.fig');
if strcmp(figPath, pngPath)
    figPath = [pngPath '.fig'];
end

try
    exportgraphics(fig, pngPath, varargin{:});
catch ME
    warning('save_fig_both: exportgraphics failed for %s: %s', pngPath, ME.message);
end

try
    savefig(fig, figPath);
catch ME
    warning('save_fig_both: savefig failed for %s: %s', figPath, ME.message);
end
end
