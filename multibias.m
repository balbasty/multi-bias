function [varargout] = multibias(x,opt)
% Estimate multiple bias fields.
%
%   This algorithm assumes a generative model of the data where multiple
%   channels represent multiple views of the same object modulated by a 
%   smooth, channel-specific multiplicative bias field (e.g., uncombined 
%   multi-coil MR images).
%   Bias fields are multiplicative and encoded by their log (ensuring
%   positivity).
%
% FORMAT [output] = multibias(x,[opt])
% x   - [Nc Ne]          Multi-view image of the same object (cell)
% b   - [Nx Ny Nz Nc]    Inferred multi-channel bias field
% m   - [Nx Ny Nz  1 Ne] Inferred mean image (after bias removal)
% opt - Structure of options with fields:
%       . itermax  - Maximum number of iterations             [16]
%       . tol      - Gain threshold for early stopping        [1E-3]
%       . lambda   - Regularisation factor    [Nc]            [1E5]
%       . sigma    - Noise standard deviation [Nc Ne]         [inferred]
%       . vs       - Voxel size                               [inferred/1]
%       . verbose  - Verbosity level: 0=quiet|1=print|2+=plot [1]
%       . map      - Memory map data (slower but saves RAM)   [false]
%       . armijo   - Damping factors for Gauss-Newton         [8 4 2 1]
%       . threads  - Number of threads; N|'matlab'            ['matlab']
%       . folder   - Output folder (empty=don't write)        ['']
%       . basename - Output basename (empty=input name)       ['']
%       . mask     - Mask of voxels to discard:               [false]
%                     false  = use all all voxels
%                     true   = estimate background mask
%                     array  = provided background mask
%       . output   - List of returned objects among:          [{'m','wb'}]
%                     'm':  mean image in mean space
%                     'wm': mean image warped to native spaces
%                     'b':  bias in mean space
%                     'wb': bias warped in native space
%                     'z':  log-bias in native space
%                     'r':  world mean-to-native rigid mappings

% -------------------------------------------------------------------------
% Input volumes
fnames = {};
dm     = [];
M      = [];
if nargin < 1, x        = [];                          end
if isempty(x), x        = spm_select(Inf, 'image');    end
if ischar(x),  x        = num2cell(x, 2);              end
if iscell(x),  fnames   = x;                           end
if iscell(x),  [x,dm,M] = files2dat(x, 'ro');          end
if isempty(x), error('No input provided');             end

Nc = size(x,1);
Ne = size(x,2);

% -------------------------------------------------------------------------
% Options
if nargin < 2, opt = struct; end
if ~isfield(opt, 'itermax'),  opt.itermax  = 16;              end
if ~isfield(opt, 'tol'),      opt.tol      = 1E-3;            end
if ~isfield(opt, 'lambda'),   opt.lambda   = 1E5;             end
if ~isfield(opt, 'sigma'),    opt.sigma    = NaN;             end
if ~isfield(opt, 'vs'),       opt.vs       = NaN;             end
if ~isfield(opt, 'verbose'),  opt.verbose  = 1;               end
if ~isfield(opt, 'mask'),     opt.mask     = 0;               end
if ~isfield(opt, 'map'),      opt.map      = 0;               end
if ~isfield(opt, 'coreg'),    opt.coreg    = true;            end
if ~isfield(opt, 'threads'),  opt.threads  = 'matlab';        end
if ~isfield(opt, 'armijo'),   opt.armijo   = [8 4 2 1];       end
if ~isfield(opt, 'folder'),   opt.folder   = '';              end
if ~isfield(opt, 'output'),   opt.output   = {'m' 'b'};       end
if ~isfield(opt, 'basename'), opt.basename = '';              end

% - Set number of threads
[nspm,nmatlab] = threads(opt.threads, 'both');

% - If no output argument, write result on disk
if nargout == 0 && isempty(opt.output)
    opt.output = spm_select(Inf, 'dir', 'Select output directory...');
end

% - Ensure right number of elements
opt.lambda = pad(opt.lambda(:), [Nc-numel(opt.lambda) 0], 'replicate', 'post');
opt.sigma  = pad(opt.sigma,     [Nc-size(opt.sigma,1) Ne-size(opt.sigma,2)], 'replicate', 'post');
if ~iscell(opt.mask), opt.mask = {opt.mask}; end
opt.mask   = pad(opt.mask,      [Nc-size(opt.mask,1) Ne-size(opt.mask,2)], 'replicate', 'post');

% - Co-register
opt.R = repmat(eye(4), [1 1 Nc Ne]);
if opt.coreg, opt.R = coreg(x, M); end
for c=1:Nc
for e=1:Ne
    M(:,:,c,e) = opt.R(:,:,c,e)\M(:,:,c,e);
end
end

% - Mean space
[dm0,opt.M0,opt.vs] = mean_space(M, dm, opt.vs);

% - Precompute mapping
opt.M = M;
for c=1:Nc
for e=1:Ne
    opt.M(:,:,c,e) = opt.M0\opt.M(:,:,c,e);
end
end

% - Load data in RAM
if ~opt.map, x = loadarray(x); end

% -------------------------------------------------------------------------
% Estimate noise standard deviation
for c=1:Nc
for e=1:Ne
if ~isfinite(opt.sigma(c,e))
    opt.sigma(c,e) = estimate_noise(x{c,e}, opt.verbose);
end
end
end

% -------------------------------------------------------------------------
% Compute background mask
for c=1:Nc
for e=1:Ne
if isscalar(opt.mask{c,e}) && opt.mask{c,e}
    warning('Not implemented yet'); opt.mask{c,e} = 0;
    % opt.mask = estimate_background(x, opt.verbose);
end
end
end

% -------------------------------------------------------------------------
% Initialisation
b  = zeros([dm0 Nc], 'single');   % Bias field
m  = update_mean(x, b, [], opt);  % Mean image
ll = [];                          % Log-likelihood

show_progress(x, b, m, ll, opt);
    
% -------------------------------------------------------------------------
% Main loop
llscl   = (prod(dm0)*Nc*Ne) * mean(1./opt.sigma.^2);
armijos = 1./opt.armijo;
for it=1:opt.itermax
    
    if opt.verbose, fprintf([repmat('-', [1 72]) '\n']); end
    
    if numel(armijos) >= it, opt.armijo = armijos(it); end
    
    [b,llx,llb] = update_bias(x, b, m, opt);
    b           = centre_bias(b, opt);
    m           = update_mean(x, b, m, opt);
   
    ll(end+1) = sum(llx(:)) + sum(llb(:));
    show_progress(x, b, m, ll, opt);
    if numel(ll) > 1, gain = -(ll(end)-ll(end-1))/llscl;
    else,             gain = Inf; end
    if opt.verbose, fprintf('Energy: %10.3g | Gain: %10.3g\n', ll(end)/llscl, gain); end
    if abs(gain) < opt.tol, break; end
end

% -------------------------------------------------------------------------
% Prepare output objects
varargout = prepare_output(b, m, x, opt);

% -------------------------------------------------------------------------
% Write on disk
if ~isempty(opt.folder), write_results(varargout, fnames, opt); end

threads(nspm,nmatlab);

% =========================================================================
%
%                                UPDATE
%
% =========================================================================

function m  = update_mean(x, b, m, opt)
% Closed-form update of the mean image

if opt.verbose, fprintf('Update mean\n'); end

Nx  = size(b,1);
Ny  = size(b,2);
Nz  = size(b,3);
Nc  = size(x,1);
Ne  = size(x,2);
dm0 = [Nx Ny Nz];

msk   = opt.mask;
M     = opt.M;
sigma = opt.sigma;

if isempty(m), m = zeros([dm0 Ne], 'single'); end

for it=1:2
for e=1:Ne
    g = single(0);
    H = single(0);
    for c=1:Nc
        if isempty(x{c,e}), continue; end
        is2  = 1./(sigma(c,e)^2);
        x1   = single(x{c,e}());
        dm1  = [size(x1) 1];
        dm1  = dm1(1:3);
        M1   = M(:,:,c,e);
        b1   = b(:,:,:,c);
        b1   = exp(pull(b1, M1, dm1));
        m1   = pull(m(:,:,:,e), M1, dm1);
        msk1 = msk{c,e};
        msk1 = msk1 | ~isfinite(x1) | x1 == 0;
        x1(msk1) = 0;
        b1(msk1) = 0;
        
        g = g + push(is2 * b1.*(b1.*m1 - x1), M1, dm0);
        H = H + push(is2 * b1.^2, M1, dm0);
    end
    m(:,:,:,e) = max(m(:,:,:,e) - max(H, eps('single')).\g, 0);
end
end

% -------------------------------------------------------------------------

function [b,llx,llb] = update_bias(x, b, m, opt)
% Gauss-Newton update of the c-th bias field

Nx  = size(b,1);
Ny  = size(b,2);
Nz  = size(b,3);
Nc  = size(x,1);
Ne  = size(x,2);
dm0 = [Nx Ny Nz];

lambda  = opt.lambda;
sigma   = opt.sigma;
armijo  = opt.armijo;
vs      = opt.vs;
verbose = opt.verbose;
msk     = opt.mask;
M       = opt.M;

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
    

    % --- Gradient and Hessian of the data term
    g = single(0);
    H = single(0);
    for e=1:Ne
        is2  = 1./(sigma(c,e)^2);
        x1   = single(x{c,e}());
        dm1  = [size(x1) 1];
        dm1  = dm1(1:3);
        M1   = M(:,:,c,e);
        b1   = exp(pull(b(:,:,:,c), M1, dm1));
        msk1 = msk{c,e};
        msk1 = msk1 | ~isfinite(x1) | x1 == 0;
        x1(msk1) = 0;
        b1(msk1) = 0;
        m1   = b1.*pull(m(:,:,:,e), M1, dm1);
        
        r1 = m1 - x1;                                   % Residuals
        llx(c,e) = 0.5 * is2 * sum(r1(:).^2, 'double'); % Log-likelihood
        g  = g + push(is2 * m1.*r1, M1, dm0);           % Gradient
        H  = H + push(is2 * m1.^2,  M1, dm0);           % Hessian
    end

    % --- Gradient of the prior term
    lam1   = lambda(c);                                    % Reg factor
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

Nx = size(b,1);
Ny = size(b,2);
Nz = size(b,3);
Nc = size(b,4);
Ne = size(m,5);
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
x1 = zeros([Nx Ny 1 Nc], 'single');
for c=1:Nc
    [x11,c11]  = push(x{c,1}, opt.M(:,:,c,1), [Nx Ny Nz]);
    x11(c11>0) = x11(c11>0)./c11(c11>0);
    x1(:,:,:,c) = x11(:,:,z);
    clear x11 c11
end
montage(x1, 'DisplayRange', [], 'ThumbnailSize', [Nx Ny]);
title('Observed images');

% --- Bias
subplot(nrow,2,4);
b1 = exp(b(:,:,z,:,1));
montage(b1, 'DisplayRange', [], 'ThumbnailSize', [Nx Ny]);
title('Bias fields');

drawnow

% -------------------------------------------------------------------------

function output = prepare_output(b, m, x, opt)

Nc = size(x,1);
Ne = size(x,2);

output = cell(1,numel(opt.output));
for i=1:numel(opt.output)
    switch opt.output{i}
        case 'm'    % Mean image in mean space
            output{i} = m;
        case 'wm'   % Mean image in native space
            output{i} = cell(size(x));
            for c=1:Nc
            for e=1:Ne
                if isempty(x{c,e}), continue; end
                dm = [size(x{c,e}) 1];
                dm = dm(1:3);
                output{i}{c,e} = pull(m, opt.M(:,:,c,e), dm);
            end
            end
        case 'b'    % Exponentiated bias field in mean space
            output{i} = exp(b);
        case 'wb'   % Exponentiated bias field in native space
            output{i} = cell(size(x));
            for c=1:Nc
            for e=1:Ne
                if isempty(x{c,e}), continue; end
                dm = [size(x{c,e}) 1];
                dm = dm(1:3);
                output{i}{c,e} = exp(pull(b(:,:,:,c,e), opt.M(:,:,c,e), dm));
            end
            end
        case 'z'    % Log bias field in mean space
            output{i} = b;
        case 'wz'   % Log bias field in native space
            output{i} = cell(size(x));
            for c=1:Nc
            for e=1:Ne
                if isempty(x{c,e}), continue; end
                dm = [size(x{c,e}) 1];
                dm = dm(1:3);
                output{i}{c,e} = pull(b(:,:,:,c,e), opt.M(:,:,c,e), dm);
            end
            end
        case 'bm'   % Modulated mean image in mean space
            output{i} = bsxfun(@times, exp(b), m);
        case 'wbm'  % Modulated mean image in native space
            output{i} = cell(size(x));
            for c=1:Nc
            for e=1:Ne
                if isempty(x{c,e}), continue; end
                dm = [size(x{c,e}) 1];
                dm = dm(1:3);
                output{i}{c,e} = pull(m(:,:,:,1,e), opt.M(:,:,c,e), dm);
                output{i}{c,e} = output{i}{c,e} .* exp(pull(b(:,:,:,c,e), opt.M(:,:,c,e), dm));
            end
            end
        case 'r'    % World mean-to-native rigid mappings
            output{i} = opt.R;
    end
end

% =========================================================================
%
%                                HELPERS
%
% =========================================================================

function write_results(output, fnames, opt)

Nc = size(opt.M,3);
Ne = size(opt.M,5);
Dc = ceil(log10(Nc+1));
De = ceil(log10(Ne+1));

% - Filenames / Orientation matrices
fnames = pad(fnames, [Nc-size(fnames,1) Ne-size(fnames,2)], '', 'post');
mat  = repmat(eye(4), [1 1 Nc Ne]);
mat0 = repmat(eye(4), [1 1 Nc Ne]);
for c=1:Nc
for e=1:Ne
    if isempty(fnames{c,e}), continue; end
    nii = nifti(fnames{c,e});
    mat(:,:,c,e)  = nii.mat;
    mat0(:,:,c,e) = nii.mat0;
end
end

% Basenames
if ~isempty(opt.basename)
    opt.basename  = ['_' opt.basename];
    basename      = cell(Nc,Ne);
    [basename{:}] = deal(opt.basename);
else
    basename      = fnames;
    for c=1:Nc
    for e=1:Ne
        if isempty(basename{c,e}), continue; end
        [~,basename{c,e}] = fileparts(basename{c,e});
        basename{c,e} = ['_' basename{c,e}];
    end
    end
end

% - Descriptions
description     = struct;
description.m   = 'Mean image (mean space)';
description.wm  = 'Mean image (native space)';
description.b   = 'Bias field (mean space)';
description.wb  = 'Bias field (native space)';
description.z   = 'Log bias field (mean space)';
description.wz  = 'Log bias field (native space)';
description.z   = 'Log bias field (mean space)';
description.wz  = 'Log bias field (native space)';
description.bm  = 'Modulated image (mean space)';
description.wbm = 'Modulated image (native space)';

% - Write files
for i=1:numel(output)
    if opt.output{i}(1) == 'w'
        % NATIVE SPACE
        for c=1:size(output{i},1)
        for e=1:size(output{i},2)
            if isempty(output{i}{c,e}), continue; end
            fname = basename{c,e};
            if size(output{i},1) > 1
                fname = sprintf([fname '-%0' num2str(Dc) 'd'], c);
            end
            if size(output{i},2) > 1
                fname = sprintf([fname '-%0' num2str(De) 'd'], e);
            end
            fname = [opt.output{i} fname '.nii'];
            nii         = nifti;
            nii.dat     = file_array(fullfile(opt.folder, fname), size(output{i}{c,e}), 'float32');
            nii.mat     = mat(:,:,c,e);
            nii.mat0    = mat0(:,:,c,e);
            nii.descrip = description.(opt.output{i});
            create(nii);
            nii.dat(:,:,:) = output{i}{c,e};
        end
        end
    elseif opt.output{i}(1) == 'r'
        % TRANSFORMATIONS
        R = opt.R;
        save(fullfile(opt.folder, ['r' basename{1,1} '.mat']), 'R');
    else
        % MEAN SPACE
        for c=1:size(output{i},4)
        for e=1:size(output{i},5)
            fname = basename{c,e};
            if size(output{i},1) > 1
                fname = sprintf([fname '-%0' num2str(Dc) 'd'], c);
            end
            if size(output{i},2) > 1
                fname = sprintf([fname '-%0' num2str(De) 'd'], e);
            end
            fname = [opt.output{i} fname '.nii'];
            dm0   = [size(output{i},1) size(output{i},2) size(output{i},3)];
            nii         = nifti;
            nii.dat     = file_array(fullfile(opt.folder, fname), dm0, 'float32');
            nii.mat     = opt.M0;
            nii.descrip = description.(opt.output{i});
            create(nii);
            nii.dat(:,:,:) = output{i}(:,:,:,c,e);
        end
        end
    end
end
% -------------------------------------------------------------------------
function [x,dm,M] = files2dat(x, permission)
% Create a virtual (memory-mapped) 5D array from a 2D cell array of
% filenames.
if nargin < 2, permission = 'rw'; end
Nc = size(x,1);
Ne = size(x,2);
dm = NaN(3,Nc,Ne);
M  = NaN(4,4,Nc,Ne);
% Convert filenames to file arrays
for c=1:Nc
for e=1:Ne
    if isempty(x{c,e}), continue; end
    if ~ischar(x{c,e}) && ~isa(x, 'nifti')
        dm1             = [size(x{c,e}) 1];
        dm1             = dm1(1:3);
        dm(:,c,e)       = dm1(:);
        M(:,:,c,e)      = eye(4);
        continue;
    end
    x{c,e}              = nifti(x{c,e});
    dm1                 = [size(x{c,e}.dat) 1];
    dm1                 = dm1(1:3);
    dm(:,c,e)           = dm1(:);
    M(:,:,c,e)          = x{c,e}.mat;
    x{c,e}              = x{c,e}.dat;
    x{c,e}.permission   = permission;
end
end
% -------------------------------------------------------------------------
function x = pull(x, y, dmo)
if ismatrix(y)
    dmi = [size(x) 1];
    dmi = dmi(1:3);
    if sum(sum((y-eye(4)).^2)) < 1E-5 && ...
       (nargin < 3 || dmo == dmi)
        return;
    end
    y  = warps_affine(dmo, y);
end
spm_diffeo('boundary', 1);
x = spm_diffeo('pull', x, y);
% -------------------------------------------------------------------------
function [x,c] = push(x, y, dmo)
if ismatrix(y)
    dmi = [size(x) 1];
    dmi = dmi(1:3);
    if sum(sum((y-eye(4)).^2)) < 1E-5 && ...
       (nargin < 3 || dmo == dmi)
        return;
    end
    y  = warps_affine(dmi, y);
end
if nargin < 3
    dmo = [size(x) 1];
    dmo = dmo(1:3);
end 
spm_diffeo('boundary', 1);
if nargout > 1
    [x,c] = spm_diffeo('push', single(x()), y, dmo);
else
    x     = spm_diffeo('push', single(x()), y, dmo);
end
% -------------------------------------------------------------------------
function psi = warps_affine(lat, mat)
id  = warps_identity(lat);
psi = reshape(reshape(id,[prod(lat) 3])*mat(1:3,1:3)' + mat(1:3,4)',[lat 3]);
if lat(3) == 1, psi(:,:,:,3) = 1; end
% -------------------------------------------------------------------------
function id = warps_identity(d)
id = zeros([d(:)' 3],'single');
[id(:,:,:,1),id(:,:,:,2),id(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
% -------------------------------------------------------------------------