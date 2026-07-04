function fit = fit_concentration_response(conc, metric, se, opts)
% FIT_CONCENTRATION_RESPONSE  Weighted log-log slope (beta) + bootstrap CI + R^2.
%
%   fit = fit_concentration_response(conc, metric, se, opts)
%
% Endpoint of the cross-algorithm reanalysis (jun9 design: simplified slope + CI + R^2,
% CCC/TOST rejected as n=8 overkill). beta = d log10(metric) / d log10(conc); beta~1 is
% linear, beta<1 sublinear (saturation/floor), beta>1 supralinear.
%
% conc   : [n] concentrations (MB/mL), > 0.
% metric : [n] response (QC-track count / amplitude), > 0.
% se     : [n] standard errors of metric for weighting (optional; [] = unweighted).
% opts   : .nBoot (default 2000), .seed (default 1; [] = unseeded), .ci (default [2.5 97.5]).
%
% Returns fit.beta, fit.beta_CI, fit.R2, fit.intercept, fit.n, fit.weighted.

    if nargin < 3, se = []; end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'nBoot'), opts.nBoot = 2000; end
    if ~isfield(opts,'seed'),  opts.seed  = 1;    end
    if ~isfield(opts,'ci'),    opts.ci    = [2.5 97.5]; end

    conc = conc(:); metric = metric(:);
    v = isfinite(conc) & isfinite(metric) & conc > 0 & metric > 0;
    if ~isempty(se), se = se(:); v = v & isfinite(se); end
    x = log10(conc(v)); y = log10(metric(v));
    n = numel(x);
    fit = struct('beta',NaN,'intercept',NaN,'R2',NaN,'beta_CI',[NaN NaN],'n',n,'weighted',~isempty(se));
    if n < 3, return; end

    % weights: propagate se into log space via the delta method, sd_log = se/(metric*ln10)
    if ~isempty(se)
        sdlog = se(v) ./ (metric(v) * log(10));
        w = 1 ./ max(sdlog, eps).^2;
    else
        w = ones(n,1);
    end

    [b0, b1] = local_wls(x, y, w);
    yhat  = b1*x + b0;
    ybar  = sum(w.*y) / sum(w);
    SSres = sum(w .* (y - yhat).^2);
    SStot = sum(w .* (y - ybar).^2);
    fit.beta = b1; fit.intercept = b0; fit.R2 = 1 - SSres / max(SStot, eps);

    % --- bootstrap CI (seeded for reproducibility) ---
    if ~isempty(opts.seed), rng(opts.seed); end
    betas = nan(opts.nBoot, 1);
    for k = 1:opts.nBoot
        idx = randi(n, n, 1);
        [~, bk] = local_wls(x(idx), y(idx), w(idx));
        betas(k) = bk;
    end
    fit.beta_CI = prctile(betas, opts.ci);
end

function [b0, b1] = local_wls(x, y, w)
% Weighted least squares for y = b1*x + b0.
    W = sum(w); Wx = sum(w.*x); Wy = sum(w.*y);
    Wxx = sum(w.*x.*x); Wxy = sum(w.*x.*y);
    d = W*Wxx - Wx^2;
    if abs(d) < eps, b1 = NaN; b0 = NaN; return; end
    b1 = (W*Wxy - Wx*Wy) / d;
    b0 = (Wy - b1*Wx) / W;
end
