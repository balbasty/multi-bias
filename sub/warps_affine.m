function psi = warps_affine(lat, mat)
id  = warps_identity(lat);
psi = reshape(reshape(id,[prod(lat) 3])*mat(1:3,1:3)' + mat(1:3,4)',[lat 3]);
if lat(3) == 1, psi(:,:,:,3) = 1; end