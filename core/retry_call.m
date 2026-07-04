function varargout = retry_call(nRetries, waitSec, fn)
%RETRY_CALL  Try fn up to nRetries times with waitSec pause between.
%
%   [a, b, ...] = retry_call(nRetries, waitSec, @() someNetworkCall(args))
%
%  Passes through nargout outputs from fn. On failure, pauses waitSec
%  seconds and retries, up to nRetries total attempts. Only rethrows
%  the error on the final attempt.
%
%  Use for VADA reads, load() on network paths, dir() on unreliable
%  network folders, etc.

for attempt = 1:nRetries
    try
        [varargout{1:nargout}] = fn();
        if attempt > 1
            fprintf('    [retry_call: succeeded on attempt %d]\n', attempt);
        end
        return;
    catch ME
        if attempt == nRetries
            fprintf('    [retry_call: exhausted %d attempts; rethrowing]\n', nRetries);
            rethrow(ME);
        end
        fprintf('    [retry_call: attempt %d/%d failed (%s); waiting %ds]\n', ...
            attempt, nRetries, ME.message, waitSec);
        pause(waitSec);
    end
end
end
