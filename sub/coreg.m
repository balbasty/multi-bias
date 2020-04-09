function R = coreg(x,M,varargin)
% Initial co-registration of input files.
%
% FORMAT M = coreg(x,M,[flags])
% x - [Nc Ne] cell of numeric|file arrays or filenames
% M - [4 4 Nc Ne] array of orientation matrices

Nc = size(x,1);
Ne = size(x,2);

% - Co-register to first volume
ref = [];
q   = NaN(6,Nc*Ne);
for i=1:numel(x)
    if isempty(x{i}), continue; end
    dat  = loadarray(x{i}, 'uint8', true, 0.9999);
    x{i} = struct('uint8', dat, 'mat', M(:,:,i));
    if isempty(ref)
        ref    = x{i};
        q(:,i) = 0;
        continue;
    end
    mov    = x{i};
    q1     = spm_coreg(ref, mov, varargin{:});
    q(:,i) = q1(:);
    x{i}   = NaN;
end

% - Find barycentric space
B = affine_basis('SE(3)');
R = NaN(4,4,Nc*Ne);
for i=1:numel(x)
    if isempty(x{i}), continue; end
    R(:,:,i) = spm_matrix(q(:,i)');
    for j=1:6 
        q(j,i) = trace(B(:,:,j)*logm(R(:,:,i))');
    end
end
R0 = spm_dexpm(mean(q,2,'omitnan'), B);

% - Apply transformation
for i=1:size(M,3)
    R(:,:,i) = R(:,:,i)/R0; % R0*(R(:,:,i)\M(:,:,i));
end
R = reshape(R, 4, 4, Nc, Ne);