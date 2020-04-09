function [x,mat,slope,inter] = loadarray(dat, otype, rnd, scale)
% Load / Convert / Scale data arrays in an (hopefully) optimal way.
%
% FORMAT [x,mat,slope,inter] = loadarray(dat, dtype, rnd, scale)
%
% dat   - fname or nifti or file_array or MATLAB array
% dtype - output data type (MATLAB style)                         [''=same]
% rnd   - if true: add noise according to data type precision       [false]
% scale - scale data to use the full dynamic range                  [false]
%         . If value in (0,1): percentile cutoff
%
% x     - MATLAB array
% mat   - orientation matrix
% slope - scaling slope
% inter - scaling intercept
%
% ----
% Note that spm_coreg uses a similar scheme to preprocess the data, with:
% . dtype = 'uint8'
% . rnd   = true
% . scale = 0.9999

if nargin < 2, otype = '';    end
if nargin < 3, rnd   = false; end
if nargin < 4, scale = false; end  

if iscell(dat)
    x     = cell(size(dat));
    mat   = cell(size(dat));
    slope = cell(size(dat));
    inter = cell(size(dat));
    for i=1:numel(dat)
        [x{i}, mat{i}, slope{i}, inter{i}] = loadarray(dat{i}, otype, rnd, scale);
    end
    return;
end

% - Dictionary of datatypes
% -------------------------------------------------------------------------
dtypes = datatypes;
mtypes = cellfun(@func2str, {dtypes.conv}, 'UniformOutput', false);
[dtypes.matlab] = deal(mtypes{:});
[dtypes(~[dtypes.supported]).matlab] = deal('unknown');

% - Read storage info
% -------------------------------------------------------------------------
% mat    <- orientation matrix   (assuming there is only one)
% itypes <- input data type(s)
% (file arrays can be formed by concatenating several other file arrays; 
%  there can therefore be more than one input type)
mat = [];
if ischar(dat), dat = nifti(dat); end
if isa(dat, 'nifti')
    mat = dat.mat;
    dat = dat.dat;
end
if isa(dat, 'file_array')
    if isempty(mat)
        fnames = {dat.fname};
        mat    = nifti(fnames);
        mat    = mat.mat;
    end
    sdat   = struct(dat);
    itypes = sdat.dtype;
    itypes = dtypes([dtypes.code] == itypes);
elseif isnumeric(dat)
    itypes = dtypes(strcmpi(class(dat), {dtypes.matlab}));
    mat    = eye(4);
else
    error('Unknown input type.');
end

% - Find output type
% -------------------------------------------------------------------------
if isempty(otype)
    if isnumeric(dat)
        otype = class(dat);
    else
        if numel(unique(itypes.code)) == 1
            otype = func2str(itypes(1).conv);
        else
            isint = true;
            if any([dat.scl_slope] ~= 1) || any([dat.scl_inter] ~= 0) || ...
               any(~(itypes.isint | itypes.code == 1))
                isint = false;
            end
            if isint
                omn = min(itypes.min);
                omx = max(itypes.max);
                if omn == 0 && omx == 1
                    otype = 'logical';
                else
                    otype = 'int';
                    if omn >= 0, otype = ['u' otype]; end
                    nbits = ceil(max(itypes.size)*8);
                    otype = [otype num2str(nbits)];
                end
            else
                if max(itypes.size) > 4, otype = 'double';
                else,                    otype = 'single'; end
            end
        end
            
    end
end
if any(strcmpi(otype, {'single' 'double' 'logical'}))
    scale = false;
end
otype = dtypes(strcmpi(otype, {dtypes.matlab}));
omn   = otype.min;
omx   = otype.max;

% - 1st pass: compute min max
% -------------------------------------------------------------------------
if scale
    imn =  Inf;
    imx = -Inf;
    if isnumeric(dat)
        imn = min(dat(:));
        imx = max(dat(:));
    else
        for i=1:numel(sdat)
            dat1 = file_array(sdat(i));
            acc  = itypes(i).isint * dat1.scl_slope;
            imn  = min(imn, min(dat1(:)));
            imx  = max(imx, max(dat1(:)) + abs(acc));
        end
    end
end

% - 2nd pass: cutoff histogram at a given percentile
% -------------------------------------------------------------------------
if scale && scale < 1
    nh = 2048;
    h  = zeros(nh,1);
    if isnumeric(dat)
        dat1 = dat();
        dat1 = dat1(isfinite(dat1));
        dat1 = round((dat1+((imx-imn)/(nh-1)-imn))*((nh-1)/(imx-imn)));
        h    = h + accumarray(dat1,1,[nh 1]);
    else
        for i=1:numel(sdat)
            dat1 = file_array(sdat(i));
            dat1 = dat1();
            dat1 = dat1(isfinite(dat1));
            dat1 = round((dat1+((imx-imn)/(nh-1)-imn))*((nh-1)/(imx-imn))); 
            h  = h + accumarray(dat1,1,[nh 1]);
        end
    end
    clear dat1
    tmp = [find(cumsum(h)/sum(h)>scale); nh];   % (spm_coreg: scale = 0.9999)
    imx  = (imn*nh-imx+tmp(1)*(imx-imn))/(nh-1);
end
 
% - Load data
% -------------------------------------------------------------------------
if rnd
    st = rand('state');
    rand('state',100);
end
if isnumeric(dat)
    x    = double(dat);
    if isinteger(dat) && rnd
        x = x + rand(size(x));
    end
    if scale
        x = (x-imn)*((omx-omn)/(imx-imn))+omn;
    end
    if otype.isint
        x = round(x);
    end
    x = max(min(x,omx),omn);
    x = otype.conv(x);
else
    x = zeros(size(dat), otype.matlab);
    for i=1:numel(sdat)
        dat1 = file_array(sdat(i));
        x1   = dat1();
        if itypes(i).isint && rnd
            x1 = x1 + dat1.scl_slope*rand(size(x1));
        end
        if scale
            x1 = (x1-imn)*((omx-omn)/(imx-imn))+omn;
        end
        if otype.isint
            x1 = round(x1);
        end
        x1 = max(min(x1,omx),omn);
        x1 = otype.conv(x1);
        S  = struct;
        S.type = '()';
        S.subs = arrayfun(@(p,d) p:(p+d-1), sdat(i).pos, sdat(i).dim, 'UniformOutput', false);
        x = subsasgn(x, S, x1);
    end
end
if rnd
    rand('state',st);
end

if scale
    slope = (imx-imn)/(omx-omn);
    inter = imn - omn * slope;
else
    slope = 1;
    inter = 0;
end