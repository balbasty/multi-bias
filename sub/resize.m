def [x, mat] = resize(x, factor, mat, anchor)
% FORMAT [x, mat] = resize(x, factor, mat, anchor)
% x      - 3D volume with optional 4th channel dimension
% factor - Resizing factor (> 1 = more/smaller voxels, < 1 = fewer/larger voxels)
% mat    - Input orientation matrix (optional)
% anchor - 'c' (centers, default) or 'e' (edges) or 'f' (first) or 'l' (last)
% 
% RETURNS
% x      - Resized volume with size floor(input_size * factor)
% mat    - Output orientation matrix
  
if nargin < 4, anchor = 'c'; end
if nargin < 3, mat    = []; end

isize = [size(x) 1];
isize = isize(1:3);
osize = floor(isize * factor);
factor = osize / isize;

g = zeros([osize 3], 'single');
switch anchor(1)
case 'c':
    scale = (isize - 1) / (osize - 1);
    shift = 1 - scale;
case 'e':
    scale = isize / osize;
    shift = 0.5 * (1 - scale);
case 'f':
    scale = isize / osize;
    shift = 1 - scale;
case 'l':
    scale = isize / osize;
    shift = isize - osize * scale;
otherwise:
    error('Unknown anchor.');
end

[g(:,:,:,1), g(:,:,:,2), g(:,:,:,3)] = meshgrid(single(1:osize(1)) * scale(1) + shift(1), ...
                                                single(1:osize(2)) * scale(2) + shift(2), ...
                                                single(1:osize(3)) * scale(3) + shift(3));

x = pull(x, g)

if nargout > 1
    scale = diag(scale);
    shift = shift(:) .* ones(3, 1);
    trf = [scale shift
           0 0 0 1];
    mat = mat * trf;
end
