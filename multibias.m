function [b,m] = multibias(x,opt)
% Estimate multiple bias fields.
%
%   This algorithm assumes a generative model of the data where multiple
%   channels represent multiple views of the same object modulated by a 
%   smooth, channel-specific multiplicative bias field (e.g., uncombined 
%   multi-coil MR images).
%   Bias fields are multiplicative and encoded by their log (ensuring
%   positivity).
%
% FORMAT [b,m] = multibias(x,[opt])
% x   - [Nx Ny Nz Nc Ne] Multi-view image of the same object
% b   - [Nx Ny Nz Nc]    Inferred multi-channel bias field
% m   - [Nx Ny Nz  1 Ne] Inferred mean image (after bias removal)
% opt - Structure of options with fields:
%       . itermax - Maximum number of iterations            [16]
%       . tol     - Gain threshold for early stopping       [1E-3]
%       . lambda  - Regularisation factor per channel       [1E3]
%       . sigma   - Noise standard deviation per channel    [inferred]
%       . vs      - Voxel size                              [inferred/1]
%       . verbose - Verbosity level 0=quiet|1=print|2=plot  [1]
%       . map     - Memory map data (slower but saves RAM)  [false]
%       . armijo  - Damping factors for Gauss-Newton        [8 4 2 1]
%       . threads - Number of threads N|'matlab'            ['matlab']
%       . folder  - Output folder.                          ['']
%       . mask    - Mask of voxels to discard:
%                   [false] = use all all voxels
%                    true   = estimate and mask out background
%                    array  = mask out using provided background mask
%
% Note that:
% . if x is a [Nc Ne] cell array of filenames, they will be concatenated to
%   for a [Nx Ny Nz Nc Ne] numeric array.
% . if x is empty, the user can select channel files through the GUI. This
%   mode does not currently allow for multiple contrasts.

% -------------------------------------------------------------------------
% Input volumes
vs     = [];
fnames = {};
if nargin < 1, x      = [];                          end
if isempty(x), x      = spm_select(Inf, 'image');    end
if ischar(x),  x      = num2cell(x, 2);              end
if iscell(x),  fnames = x;                           end
if iscell(x),  [x,vs] = files2dat(x, 'ro');          end
if isempty(x), error('No input provided');           end

Nx  = size(x,1);
Ny  = size(x,2);
Nz  = size(x,3);
Nc  = size(x,4);
Ne  = size(x,5);
lat = [Nx Ny Nz];

% -------------------------------------------------------------------------
% Options
if nargin < 2, opt = struct; end
if ~isfield(opt, 'itermax'), opt.itermax = 16;              end
if ~isfield(opt, 'tol'),     opt.tol     = 1E-3;            end
if ~isfield(opt, 'lambda'),  opt.lambda  = 1E5;             end
if ~isfield(opt, 'sigma'),   opt.sigma   = NaN;             end
if ~isfield(opt, 'vs'),      opt.vs      = NaN;             end
if ~isfield(opt, 'verbose'), opt.verbose = 1;               end
if ~isfield(opt, 'mask'),    opt.mask    = 0;               end
if ~isfield(opt, 'map'),     opt.map     = 0;               end
if ~isfield(opt, 'threads'), opt.threads = 'matlab';        end
if ~isfield(opt, 'armijo'),  opt.armijo  = [8 4 2 1];       end
if ~isfield(opt, 'output'),  opt.output  = '';              end

% - Set number of threads
[nspm,nmatlab] = threads(opt.threads, 'both');

% - If no output argument, write result on disk
if nargout == 0 && isempty(opt.output)
    opt.output = spm_select(Inf, 'dir', 'Select output directory...');
end

% - Set inferred voxel size
if isscalar(opt.vs) && ~isfinite(opt.vs)
    if isempty(vs), opt.vs = 1;
    else,           opt.vs = vs; end
end

% - Ensure right number of elements
opt.lambda = pad(opt.lambda(:), [Nc - numel(opt.lambda) 0], 'replicate', 'post');
opt.sigma  = pad(opt.sigma(:),  [Nc - numel(opt.sigma)  0], 'replicate', 'post');
opt.vs     = pad(opt.vs(:)',    [3  - numel(opt.vs)     0], 'replicate', 'post');

% - Load data in RAM
if ~opt.map, x = loadsametype(x); end

% -------------------------------------------------------------------------
% Estimate noise standard deviation
if any(~isfinite(opt.sigma))
    % 1) Estimate noise sd in each channel and contrast
    % 2) Weighted geometric mean across contrasts
    [sigma,mu] = estimate_noise(x, opt.verbose);
    sigma      = exp(sum(log(sigma).*mu,2)./sum(mu,2));
    opt.sigma(~isfinite(opt.sigma)) = sigma(~isfinite(opt.sigma));
end

% -------------------------------------------------------------------------
% Compute background mask
if isscalar(opt.mask) && opt.mask
    warning('Not implemented yet'); opt.mask = 0;
    % opt.mask = estimate_background(x, opt.verbose);
end
% Make it a mask of kept voxels (instead of discarded voxels)
opt.mask = ~opt.mask;

% -------------------------------------------------------------------------
% Initialisation
b  = zeros([lat Nc], 'single');   % Bias field
m  = update_mean(x, b, opt);      % Mean image
ll = [];                          % Log-likelihood

show_progress(x, b, m, ll, opt);
    
% -------------------------------------------------------------------------
% Main loop
llscl   = (Nx*Ny*Nz*Nc*Ne) * mean(1./opt.sigma.^2);
armijos = 1./opt.armijo;
for it=1:opt.itermax
    
    if opt.verbose, fprintf([repmat('-', [1 72]) '\n']); end
    
    if numel(armijos) >= it, opt.armijo = armijos(it); end
    
    [b,llx,llb] = update_bias(x, b, m, opt);
    b           = centre_bias(b, opt);
    m           = update_mean(x, b, opt);
   
    ll(end+1) = sum(llx(:)) + sum(llb(:));
    show_progress(x, b, m, ll, opt);
    if numel(ll) > 1, gain = -(ll(end)-ll(end-1))/llscl;
    else,             gain = Inf; end
    if opt.verbose, fprintf('Energy: %10.3g | Gain: %10.3g\n', ll(end)/llscl, gain); end
    if abs(gain) < opt.tol, break; end
end

% -------------------------------------------------------------------------
% Exponentiate bias field
b = exp(b);

% -------------------------------------------------------------------------
% Write on disk
if ~isempty(opt.output), write_results(b, m, fnames, opt); end

threads(nspm,nmatlab);

% =========================================================================
%
%                                UPDATE
%
% =========================================================================

function m  = update_mean(x, b, opt)
% Closed-form update of the mean image

if opt.verbose, fprintf('Update mean\n'); end

Nx  = size(x,1);
Ny  = size(x,2);
Nz  = size(x,3);
Nc  = size(x,4);
Ne  = size(x,5);

m = zeros([Nx Ny Nz 1 Ne], 'single');
z = zeros([Nx Ny Nz],      'single');
for c=1:Nc
    is2 = 1./(opt.sigma(c)^2);
    b1  = exp(b(:,:,:,c));
    z   = z + is2 * b1.^2;
    for e=1:Ne
        m(:,:,:,1,e) = m(:,:,:,1,e) + is2 * b1 .* single(x(:,:,:,c,e));
    end
end
m = bsxfun(@rdivide, m, z);

% -------------------------------------------------------------------------

function [b,llx,llb] = update_bias(x, b, m, opt)
% Gauss-Newton update of the c-th bias field

Nc  = size(x,4);
Ne  = size(x,5);

lambda  = opt.lambda;
sigma   = opt.sigma;
armijo  = opt.armijo;
vs      = opt.vs;
verbose = opt.verbose;
mask    = opt.mask;

if verbose
    ndigitsc = ceil(log10(Nc+1));
    fprintf('Update bias');
    if Nc > 1
        fprintf(': [ ');
        fprintf(['%0' num2str(ndigitsc) 'd '], 0);
        fprintf(']');
    end
end

llx = zeros(Nc,Ne);
llb = zeros(Nc,1);
    
for c=1:Nc

    if verbose && Nc > 1
        fprintf('\b');
        fprintf(repmat('\b', [1 ndigitsc + 1]));
        fprintf(['%0' num2str(ndigitsc) 'd '], c);
        fprintf(']');
    end
    
    lam1 = lambda(c);        % Regularisation factor
    is2  = 1./(sigma(c)^2);  % Inverse variance
    b1   = exp(b(:,:,:,c));  % Exponentiated bias field

    % --- Gradient and Hessian of the data term
    g = single(0);
    H = single(0);
    for e=1:Ne
        m1 = b1 .* m(:,:,:,1,e) .* mask;                % Modulated mean
        r1 = m1 - single(x(:,:,:,c,e)) .* mask;         % Residuals
        llx(c,e) = 0.5 * is2 * sum(r1(:).^2, 'double'); % Log-likelihood
        g  = g + is2 * m1 .* r1;                        % Gradient
        H  = H + is2 * m1 .^ 2;                         % Hessian
    end

    % --- Gradient of the prior term
    b1     = b(:,:,:,c);                                   % Log bias field
    mom    = spm_field('vel2mom', b1, [vs 0 0 lam1]);      % Momentum
    g      = g + mom;                                      % Gradient
    llb(c) = 0.5 * sum(mom(:).*b1(:), 'double');           % Log-likelihood
    
    % --- Update step
    b1  = b1 - armijo * spm_field(H, g, [vs 0 0 lam1 2 2]);
    b(:,:,:,c) = b1;
    
end
if verbose, fprintf('\n'); end

% -------------------------------------------------------------------------

function b = centre_bias(b, opt)
% Zero-centre the bias fields

if opt.verbose, fprintf('Centre bias\n'); end

sumb = single(0);
for c=1:size(b,4)
    sumb = sumb + opt.lambda(c) * b(:,:,:,c);
end
sumb = sumb / sum(opt.lambda);
for c=1:size(b,4)
    b(:,:,:,c) = b(:,:,:,c) - sumb;
end

% -------------------------------------------------------------------------

function show_progress(x, b, m, ll, opt)

if opt.verbose < 2, return; end
   
figname = 'MultiBias';
f = findobj('Type', 'Figure', 'Name', figname);
if isempty(f)
    f = figure('Name', figname, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);   
clf(f);

if opt.verbose > 2, nrow = 2;
else,               nrow = 1; end

Nx = size(x,1);
Ny = size(x,2);
Nz = size(x,3);
Nc = size(x,4);
Ne = size(x,5);
z  = ceil(Nz/2);

% --- Mean image
subplot(nrow,2,1);
m1 = reshape(m(:,:,z,1,:), [Nx Ny 1 Ne]);
mu = reshape(mean(reshape(m1, [], Ne)), [1 1 1 Ne]);
m1 = bsxfun(@rdivide, m1, mu);
montage(m1, 'DisplayRange', [], 'ThumbnailSize', [Nx Ny]);
colormap('gray');
title('Mean image');

% --- Log-likelihood
subplot(nrow,2,2);
plot(-ll);
title('Log-likelihood');

if opt.verbose == 2, drawnow; return; end

% --- Data
subplot(nrow,2,3);
x1 = reshape(x(:,:,z,:,1), [Nx Ny 1 Nc]);
montage(x1, 'DisplayRange', [], 'ThumbnailSize', [Nx Ny]);
title('Observed images');

% --- Bias
subplot(nrow,2,4);
b1 = reshape(exp(b(:,:,z,:,1)), [Nx Ny 1 Nc]);
montage(b1, 'DisplayRange', [], 'ThumbnailSize', [Nx Ny]);
title('Bias fields');

drawnow

% =========================================================================
%
%                                HELPERS
%
% =========================================================================

function write_results(b, m, fnames, opt)

Nx = size(b,1);
Ny = size(b,2);
Nz = size(b,3);
Nc = size(b,4);
Ne = size(m, 5);

ndigitsc = ceil(log10(Nc+1));
ndigitse = ceil(log10(Ne+1));

% - Create filenames for bias fields
bnames   = cell(1,Nc);  % Filenames
bmats    = cell(1,Nc);  % Orientation matrices
bdescrip = cell(1,Nc);  % Description
if size(fnames,1) == Nc
    for c=1:Nc
        nii           = nifti(fnames{c,1});
        [~, fname, ~] = fileparts(fnames{c,1});
        bnames{c}     = ['b' fname '.nii'];
        bmats{c}      = nii.mat;
        bdescrip{c}   = nii.descrip;
    end
elseif ~isempty(fnames)
    nii           = nifti(fnames{1});
    [~, fname, ~] = fileparts(fnames{1});
    for c=1:Nc
        bnames{c}   = sprintf(['b' fname '-%0' num2str(ndigitsc) 'd.nii'], c);
        bmats{c}    = nii.mat;
        bdescrip{c} = nii.descrip;
    end
else
    for c=1:Nc
        bnames{c}   = sprintf(['b-%0' num2str(ndigitsc) 'd.nii'], c);
        bmats{c}    = eye(4);
        bdescrip{c} = sprintf('Bias field [%d]', c);
    end
end

% - Create filenames for mean images
mnames   = cell(1,Ne);  % Filenames
mmats    = cell(1,Ne);  % Orientation matrices
mdescrip = cell(1,Ne);  % Description
if size(fnames,2) == Ne
    for e=1:Ne
        nii           = nifti(fnames{1,e});
        [~, fname, ~] = fileparts(fnames{1,e});
        mnames{e}     = ['m' fname '.nii'];
        mmats{e}      = nii.mat;
        mdescrip{e}   = nii.descrip;
    end
elseif ~isempty(fnames)
    nii           = nifti(fnames{1});
    [~, fname, ~] = fileparts(fnames{1});
    for e=1:Ne
        mnames{e}   = sprintf(['m' fname '-%0' num2str(ndigitse) 'd.nii'], e);
        bmats{e}    = nii.mat;
        bdescrip{e} = nii.descrip;
    end
else
    for e=1:Ne
        mnames{e} = sprintf(['m-%0' num2str(ndigitse) 'd.nii'], e);
        mmats{e}    = eye(4);
        mdescrip{e} = sprintf('Mean image [%d]', e);
    end
end

% - write files
nii0     = nifti;
nii0.dat = file_array(fullfile(opt.output, bnames{1}), [Nx Ny Nz], 'float32');
for c=1:Nc
    nii           = nii0;
    nii.dat.fname = fullfile(opt.output, bnames{c});
    nii.mat       = bmats{c};
    nii.descrip   = bdescrip{c};
    create(nii);
    nii.dat(:,:,:) = b(:,:,:,c);
end
for e=1:Ne
    nii           = nii0;
    nii.dat.fname = fullfile(opt.output, mnames{e});
    nii.mat       = mmats{e};
    nii.descrip   = mdescrip{e};
    create(nii);
    nii.dat(:,:,:) = m(:,:,:,1,e);
end

function x = loadsametype(x)
% Load file array while retaining on-disk datatypes
% (by default, file arrays are loaded as 'double')

if ~isa(x, 'file_array'), return; end

xt  = unique([struct(x).dtype]);
dt  = datatypes;
xt  = dt(ismember([dt.code], xt));
if any(~[xt.isint])
    xt = xt(~[xt.isint]);
    converter = xt([xt.size] == max([xt.size])).conv;
elseif any(~[xt.unsigned])
    converter = dt(~[dt.unsigned] & [dt.int] & [dt.size] == max([xt.size])).conv;
else
    converter = xt([xt.size] == max([xt.size])).conv;
end
x = converter(x());

function [xo,vs] = files2dat(x, permission)
% Create a virtual (memory-mapped) 5D array from a 2D cell array of
% filenames.
if nargin < 2, permission = 'rw'; end
Nc = size(x,1);
Ne = size(x,2);
vs = [];
% Convert filenames to file arrays
for c=1:Nc
for e=1:Ne
    x{c,e} = nifti(deblank(x{c,e}));
    if isempty(vs), vs = sqrt(sum(x{c,e}.mat(1:3,1:3).^2)); end
    x{c,e} = x{c,e}.dat;
    x{c,e}.permission = permission;
end
end
% Concatenate file arrays
xo = cell(1,Nc);
for c=1:Nc
    xo{c} = cat(5, x{c,:});
end
xo = cat(4, xo{:});