function x = pull(x, y, dmo)
if ismatrix(y)
    dmi = [size(x) 1];
    dmi = dmi(1:3);
    if sum(sum((y-eye(4)).^2)) < 1E-5 && ...
       (nargin < 3 || all(dmo == dmi))
        return;
    end
    y  = warps_affine(dmo, y);
end
spm_diffeo('boundary', 1);
x = spm_diffeo('pull', x, y);