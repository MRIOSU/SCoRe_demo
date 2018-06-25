clear;

%% Add relevant folders to the path
restoredefaultpath
if ( ~isempty( strfind(pwd,'\')) ) % Windows('contains' doesn't work in Linux )
addpath('.\bFISTA');
addpath('.\UWT');
addpath('.\coilSensitivity');
else
addpath('./bFISTA'); % Linux
addpath('./UWT'); 
addpath('./coilSensitivity');
end


%% Data parameters
p.nStd0  = 1; % Noise power after noise prewhiting
p.yMx    = 50; % Normalize the data to this max value
p.isoRes  = 1; % zero-pad PE to make the spatial resolution isotropic

%% sensitivity estimation parameters
p.fil    = [6,6,6];  % size of kernal for eSPIRiT
p.ACSsz  = [64,64,30]; % size of the k-space block used for eSPIRiT
p.eSRT   = 0.95;
p.avgPhs = 1; % Assign the phase to time-average image. 1: yes, 0: no

%% Sparsifying transform parameters
p.N        = 1;  % level of decomposition
p.dict = {{'db1','db1','db1'}}; % 'fd'; % concatenation of n wavelet transforms {'db1','db2',....'dbn'},                .
p.bGrp = ones(1, 8); % band grouping for SCoRe, each group is threholded similarly


%% SCoRe parameters
p.sLP      = 1; % the distribution may be skewed along real or imaginary axis;
 
p.nFtr     = 1; % Global scaling of the noise power
p.sFtri    = 2; % Global "initial" scaling of the composite sparsity term
p.sFtrf    = 1; % Global "final" scaling of the composite sparsity term

p.nTol	   = [5e-3, 1/10]; % Noise toleratance; lambda = 1/(L1 + noise tolerance)
p.lRes     = [0.3, 5]; % Restriction on lambda values; for p.lRes(1) interations, keep max(lambda)<= p.lRes(2)*meanLambda


%% bFISTA paramters
p.oIter = 10;  % Total outer iteration
p.iIter    = 8; % Inner FISTA iterations
p.minIter  = 1/2; % Minimum number of inner iterations (as a fraction of p.iIter) to run for each outer iteration
p.fstIter  = 2; % for oIter = 1, the iIter = p.iIter x p.fstIter
p.stpThrsh = 2e-6; % Stopping threshold
p.L1       = []; % Lipschitz constant for the fidelity term; will be calculated if left empty



%% Load data and perform recon
filename = 'input1_4CH.mat';
load(['.\data\' filename]);%dims of data is fixed as [RO E1 E2 CHA SLC PHS other], Noise power = 1 for each channel
[y,samp] = dataAjst(data,p,param);
p.param = param;

%%%%%%%%%%%%%%%%%       Coil sensitivity estimation     %%%%%%%%%%%%%%%%%%%
tic;[S,x0]  = coilSen(y, samp, p);toc;

%%%%%%%%%%%%%%%%%%%%%%%%%  Data  Scaling   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crop data and find matrix size
[y, yMat, samp, sampInd, sclFctr, nSize] = dataScal(y, p);
p.n = [nSize(1:end-2),nSize(end)]; % Image size (remove the coil_dim in k-space)
p.R = numel(samp(floor(end/2)+1,:,:))/sum(sum(samp(floor(end/2)+1,:,:))); % Net acceleration
dim_image = numel(p.n); % image dimmension

%%%%%%%%%%%%%%%%%%%%%%%%%  Initialization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.nStd  = p.nStd0*sclFctr; % noise power after scaling
p.snr   = 10*log10(sum(abs(y).^2)/(numel(y)*p.nStd^2)); % SNR in db
p.sFtri = max(0.5, exp(p.snr/10-3));

%Forward operator and sparsifying transform
p.A  = @(x) funA (x, sampInd, S, nSize);
p.At = @(x) funAt(x, sampInd, S, nSize);
[p.U, p.Ut, p.bInd] = wavWrap(randn(p.n), p);

%  Estimate noise in  image domain
tmp1 = p.At(p.nStd/sqrt(2)*randn(size(y)) + 1j*p.nStd/sqrt(2)*randn(size(y)));
p.nStdI = p.R*std(tmp1(:)); % Noise in the image domain; p.R for amplification.
clear tmp1;

% Initial x0 and lm0
x0  = x0*sclFctr;
x0  = repmat(x0,[1,1,p.n(end)]);%Insert frame dim
x0  = x0 + 0.010*max(abs(x0(:)))*randn(size(x0)); % inital image guess

%x0 = p.At(y);

p.lmb0  = 0.25/mean(abs(p.U(x0)));   % Initial /lambda guess

%%%%%%%%%%%%%%%%%%%%%%%%%  bFISTA slover   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;[tmp,lmb,L1] = bFISTA(y, x0, p);toc % reconstrution with data weighting
%
%% reconstructed image xHat
xHat = 1/sclFctr * tmp;
%display the image
figure;
for n = 1:2
    for fr = 1:size(xHat,3)
        imagesc(abs(xHat(:,:,fr)),[0,8*mean(abs(xHat(:)))]); colormap(gray); axis image; axis off; pause(0.001*p.param.TRes);
    end
end
p.lmb = lmb(p.bInd(1:end-1)+1);
pOut = p; pOut.A = []; pOut.At =[]; pOut.U = []; pOut.Ut =[];
save(['./data/output' filename(6:end)],'xHat','pOut');
close all;
