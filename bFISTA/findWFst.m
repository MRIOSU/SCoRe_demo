function [wt, L] = findWFst(coefX, param)
% Estimate tuning parameters

clear Ux;
bInd   = param.bInd;
nCoef  = bInd(end);

nStd   = param.nStd;
nStdW  = param.nStdI; % Noise in image/transform domain


Nlim   = param.Nlim;
wt     = zeros(nCoef,1);

avCoefX= max(mean(abs(coefX(:)))*param.nTol(1), nStdW*param.nTol(2));

sLP    = param.sLP;
bGrp   = param.bGrp;
bGrp   = [0, cumsum(bGrp)];

nFtr   = 1 * param.nFtr * nStd^2;
sFtr   = 2 * param.sFtr;

    % ====SCoRe=================================
    L   = ones(numel(bGrp)-1, 1);
    for j = 1:(numel(bGrp)-1)
        tmp = coefX(bInd(bGrp(j)+1)+1:bInd(bGrp(j+1)+1));
        if  sLP == 1
            LTmp = [1/(sum(abs(real(tmp)))/numel(tmp) + avCoefX), 1/(sum(abs(imag(tmp)))/numel(tmp) + avCoefX)];
            L(j) = 1/2*sqrt((1+(min(LTmp(1),LTmp(2))/max(LTmp(1),LTmp(2)))^2))*min(LTmp(1),LTmp(2)); % Find a compromise between circular laplacian and skewed laplacian
            wt(bInd(bGrp(j)+1)+1:bInd(bGrp(j+1)+1)) = nFtr * sFtr *L(j);
        else
            L(j) = 1/(sum(abs(tmp))/numel(tmp) + avCoefX); % The numerator becomes "2" for complex Laplacian prior; 1e-3 some random scaling to keep numbers small
            wt(bInd(bGrp(j)+1)+1:bInd(bGrp(j+1)+1)) = nFtr * sFtr *L(j);
        end
    end
    gmL = 1/(sum(abs(coefX(:)))/numel(coefX(:))); % "average" L
    gmwt = nFtr * sFtr *gmL;
    wt(wt < 1/Nlim*gmwt) = 1/Nlim*gmwt;
    wt(wt > Nlim*gmwt)   = Nlim*gmwt;
    L(L < 1/Nlim*gmL) = 1/Nlim*gmL;
    L(L > Nlim*gmL)   = Nlim*gmL;
    disp(['Sub-band lambda values: ' num2str(L', 4)]);

