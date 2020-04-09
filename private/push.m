function [x,c] = push(x, y, dmo)
if ismatrix(y)
    dmi = [size(x) 1];
    dmi = dmi(1:3);
    if sum(sum((y-eye(4)).^2)) < 1E-5 && ...
       (nargin < 3 || all(dmo == dmi))
        c = ones(dmi, 'single');
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