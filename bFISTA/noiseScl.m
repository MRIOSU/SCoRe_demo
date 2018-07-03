function[x] = noiseScl(x)
warning('Noise variance estimated from the k-space fringes may not be accurate, especially when SNR is very high');

y = squeeze(x); % read out, phase encode, coil, frame

% discard the central part of k-space
tmp = y([1:min(16,round(size(y,1)/8)), end-min(16,round(size(y,1)/8))+1:end],:,:,:);
tmp = tmp(:,[1:round(size(y,2)/3), end-round(size(y,2)/3)+1:end],:,:);

cStd = zeros(size(y,3),1);
 for i = 1:size(y,3)
        tmp2 = tmp(:,:,i,:);
        tmp2 = std(tmp2(tmp2~=0)); 
        cStd(i) = tmp2; % Noise in the ith coil
 end
 
nStd = prctile(sort(cStd), 10)/1.2; % Guess the noise variance from noisy coil
% nStd = median(abs(cStd))/sqrt(0.675)/2;
x= x/nStd;
