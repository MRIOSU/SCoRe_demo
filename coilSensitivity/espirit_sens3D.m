%% ESPIRiT Maps Demo
% This is a demo on how to generate ESPIRiT maps. It is based on the paper
% Uecker et. al, MRM 2013 DOI 10.1002/mrm.24751. ESPIRiT is a method that
% finds the subspace of multi-coil data from a calibration region in
% k-space using a series of eigen-value decompositions in k-space and image
% space. 
function [cI, maps] = espirit_sens3D(DATA, samp, ncalib, param)
% eigThresh_1 = param.eSRT(1);
eigThresh_2 = param.eSRT;
ksize = param.fil;
ksize = min(ksize,ncalib);

[sx,sy,sz,Nc] = size(DATA);


% crop a calibration area
calib = crop(DATA,[ncalib,Nc]);
samp = crop(samp,ncalib(1:ndims(samp)));
%%
% Display coil images: 
im = DATA;
for n = 1:3
im = ifftc(im,n);
end


%% Compute ESPIRiT EigenVectors
% Here we perform calibration in k-space followed by an eigen-decomposition
% in image space to produce the EigenMaps. 


% compute Calibration matrix, perform 1st SVD and convert singular vectors
% into k-space kernels

[k,S,dim_A] = dat2Kernel3D(calib,samp,ksize);

% estimate the noise variance
zS = max(0,(dim_A(2)-numel(S)));
S_tmp = padarray(S,[zS 0],0,'post');

seed = 100;
nStd0 = 1;
rng(100*seed); n_real = randn(dim_A);
rng(200*seed); n_im = randn(dim_A);
noise  = nStd0*complex(n_real,n_im)/sqrt(2); % noise
[~,S_noise,~] = svd(noise,'econ');
S_noise = diag(S_noise);
sigma = sqrt((mean(S(ceil(3*end/4):end))/mean(S_noise(ceil(3*end/4):end))));
clear noise;

% soft threshold lambda
lambda = S(end-1); %inital guess of lambda
IS_REAL = 0;
MSE = sure_svt(lambda, sigma, S, dim_A, IS_REAL);
% MSEcc = zeros(numel(S),1);
for jdx = 1:(numel(S)-1) % figure out the optimal lambda
    testlambda = S(end-jdx);
    testMSE = sure_svt(testlambda, (sigma), S, dim_A, IS_REAL);
%     MSEcc(end-jdx) = testMSE; 
    if testMSE < MSE
        lambda = testlambda;
        MSE = testMSE;
    end   
end
% weight = max(0,S-lambda)./(abs(S-lambda)+eps);
weight = max(0,S-lambda)./(S+eps);
idx = max(find(weight > 0));


%% 
% This shows that the calibration matrix has a null space as shown in the
% paper. 
disp([num2str(100*double(idx)/numel(S(:))) '% of singular vectors are included.']);
if ( double(idx)/numel(S(:)) >= 0.7 || double(idx)/numel(S(:)) <= 0.2)
     warning('ESPiRiT parameter estimated by SURE may not be resonable!');
     kdisp = reshape(k,[ksize(1)*ksize(2)*ksize(3)*Nc,ksize(1)*ksize(2)*ksize(3)*Nc]);
     figure, subplot(211), plot([1:ksize(1)*ksize(2)*ksize(3)*Nc],S_tmp,'LineWidth',2);
     hold on, 
     plot([1:ksize(1)*ksize(2)*ksize(3)*Nc],S_tmp(idx),'r-','LineWidth',2);
     plot([idx,idx],[0,S_tmp(1)],'g--','LineWidth',2)
     legend('signular vector value','threshold')
     title('Singular Vectors')
     subplot(212), imagesc(abs(kdisp)), colormap(gray(256));
     xlabel('Singular value #');
     title('Singular vectors')
end


%%
% crop kernels and compute eigen-value decomposition in image space to get
% maps
[M,W] = kernelEig3D(k(:,:,:,:,1:idx),[sx,sy,sz]);


%%
% project onto the eigenvectors. This shows that all the signal energy
% lives in the subspace spanned by the eigenvectors with eigenvalue 1.
% These look like sensitivity maps. 
P = sum(repmat(im,[1,1,1,1,Nc]).*conj(M),4);


%%
% crop sensitivity maps 
maps = M(:,:,:,:,end);
cI = permute(P(:,:,:,1,end),[1,2,3,5,4]); % Coil combined image,(kx,ky,kz,set)
% Weight the eigenvectors with soft-senses eigen-values
weights = W(:,:,:,end) ;
weights = (weights - eigThresh_2)./(1-eigThresh_2).* (W(:,:,:,end) > eigThresh_2);
weights = -cos(pi*weights)/2 + 1/2;
maps = maps.*sqrt(repmat(permute(weights,[1 2 3 5 4]),[1 1 1 Nc 1]));



