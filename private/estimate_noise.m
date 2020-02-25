function [sigma,mu,info] = estimate_noise(x, verbose)
% Estimate the noise precision of each channel by fitting a Rice mixture.
%
% FORMAT [sigma,mu,info] = estimate_noise(x, [verbose])
%
% x       - [Nx Ny Nz Nc Ne] - Observed volume
% verbose                    - 0=quiet | [1=print]
% sigma   - [Nc Ne]          - Diagonal noise precision matrix
% mu      - [Nc Ne]          - Mean tissue signal
% info    - [Ne Ne]          - Can be used for plotting.
%
% `info` can be used for plotting the fit as:
% >> plot(info.x(:),info.p,'--',info.x(:), ...
%         info.h/sum(info.h)/info.md,'b.', ...
%         info.x(:),info.sp,'r');
%
% If the noise s.d. is known to be share across echoes, a good estimate can
% be obtained by taking the geometric mean of sigma across echoes:
% >> exp(mean(log(sigma),2))

if nargin < 2, verbose = 1; end

Nc    = size(x,4);
Ne    = size(x,5);
sigma = zeros([Nc Ne]);
mu    = zeros([Nc Ne]);
info  = struct('x',[],'h',[],'p',[],'sp',[],'md',[]);

% -------------------------------------------------------------------------
% Verbosity
if verbose
    ndigitsc = ceil(log10(Nc+1));
    ndigitse = ceil(log10(Ne+1));
    fprintf('Estimate noise');
    if Nc > 1 || Ne > 1
        fprintf(': [ ');
        if Nc > 1, fprintf(['%0' num2str(ndigitsc) 'd '], 0); end
        if Ne > 1, fprintf(['%0' num2str(ndigitse) 'd '], 0); end
        fprintf(']');
    end
end

% -------------------------------------------------------------------------
% Loop over 3D volumes
for ne=1:Ne
for nc=1:Nc
    
    if verbose
        if Nc > 1 || Ne > 1, fprintf('\b'); end
        if Nc > 1, fprintf(repmat('\b', [1 ndigitsc + 1])); end
        if Ne > 1, fprintf(repmat('\b', [1 ndigitse + 1])); end
        if Nc > 1, fprintf(['%0' num2str(ndigitsc) 'd '], nc); end
        if Ne > 1, fprintf(['%0' num2str(ndigitse) 'd '], ne); end
        if Nc > 1 || Ne > 1, fprintf(']'); end
    end
    
    
    % ---------------------------------------------------------------------
    % Guess data type
    if isa(x, 'file_array')
        % Try to find which file_array corresponds to the current volume
        dat   = struct(x);
        found = false;
        for k=1:numel(dat)
            pos = [dat(k).pos 1 1 1];
            dim = [dat(k).dim 1 1 1];
            if nc >= pos(4) && nc < pos(4) + dim(4) && ...
               ne >= pos(5) && ne < pos(5) + dim(5)
                found = true;
                break
            end
        end
        if ~found
            warning('This should not happen... Using the first array instead. ' );
            k = 1;
        end
            
        dtype = {x.dtype};
        dtype = dtype{k};
        isint = spm_type(dtype(1:(end-3)),'intt');
        slope = {x.scl_slope};
        slope = slope{k};
    else
        isint = isinteger(x);
        slope = 1;
    end
    
    % ---------------------------------------------------------------------
    % Fit Rician mixture
    xn = single(x(:,:,:,nc,ne));
    [sigma(nc,ne), mu(nc,ne), info(nc,ne)] = estimate_noise_1(xn, 2, isint, slope, nargout > 1);
    
end
end
if verbose, fprintf('\n'); end

% =========================================================================
% Helper: estimate noise from one 3D scan
function [noise, mu_val, info] = estimate_noise_1(f, K, inttype, slope, getmu)

if inttype
    f(f==max(f(:))) = 0;
    x      = 0:slope:max(f(:));
    [h,x]  = hist(f(f~=0),x);
else
    x      = (0:1023)*(max(f(:))/1023);
    f(f==max(f(:))) = 0;
    [h,x]  = hist(f(f~=0 & isfinite(f)),x);
end
[mg,nu,sd,info] = spm_rice_mixture(double(h(:)),double(x(:)),K);
noise           = min(sd);

mu_val = NaN;
if getmu && nargout>=2
    x          = -nu.^2./(2*sd.^2);
    msk        = x>-20;
    Laguerre   = exp(x(msk)/2).*((1-x(msk)).*besseli(0,-x(msk)/2)-x(msk).*besseli(1,-x(msk)/2));
    Ey         = zeros(size(sd));
    Ey( msk)   = sqrt(pi*sd(msk).^2/2).*Laguerre;
    Ey(~msk)   = nu(~msk);
    mu_val     = max(Ey);
end