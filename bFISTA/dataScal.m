function [y, yMat, samp, sampInd, sclFctr, nSize] = dataScal(yMat, p)
% input: [RO E1 E2 CHA SLC PHA others]
% output: [RO E1 E2 CHA PHA others]

% scale data 
sclFctr = p.yMx/max(abs(yMat(:)));
yMat    = yMat*sclFctr;

%    Vectorize data
samp = logical(abs(yMat)); %(kx,ky,kz,coil,frame)
sampInd = find(samp~=0);
y = yMat(sampInd);

% pesudo 3D to 2D
samp = squeeze(samp(:,:,:,1,1,:,:));
yMat = squeeze(yMat);
nSize = size(yMat);

end
