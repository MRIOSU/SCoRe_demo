function [y] =funAt(x,sampInd, S, nSize)
% works for 2D,3D,single and multiple sensitivity maps
% k-space to image
y = zeros(nSize);

%2D-kspace(maps) to pesudo 3D-kspace(maps)
if numel(nSize) == 4 %kx,ky,coil,frame
    y = permute(y,[1,2,5,3,4]); %kx,ky,kz,coil,frame
    S = permute(S,[1,2,5,3,4]); %kx,ky,kz,coil,set    
end

y(sampInd) = x;
for n = 1:3
  y = ifftc(y,n);
end
y = squeeze(1 * sum(bsxfun(@times, permute(conj(S),[1 2 3 4 6 5]), y), 4)); %kx,ky,kz,coil,frame,set; sum coil dim
