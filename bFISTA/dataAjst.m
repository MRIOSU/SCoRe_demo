function [y, samp] = dataAjst(Y, p, param)
% input: [RO E1 E2 CHA SLC PHS others]
% output: [RO E1 E2 CHA SLC PHS others] (dim_slc = 1)

y = Y;

%make PE even(first PE direction)
if rem(size(y,2), 2) == 1
    y(:,end+1,:,:,:,:,:)=0;    % only use this to make PE even
end

if p.isoRes == 1 % make resolution isotropic
    zpy = max(round(1/2 * (param.FOV(2)/param.FOV(1) * size(y,1) - size(y,2))), 0);
    zpz = 0;
    y = padarray(y, [0, zpy, zpz], 0, 'both');
end        

samp = logical(abs(y));
samp = samp(:,:,:,1,1,:,:); %set dim_CHA, dim_SLC = 1

