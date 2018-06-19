function [x, lmb, G, L1] = bFISTA(b, x0, param)
%-----------------------------------------------------------------------
% Adhmad.46@osu.edu
% Last modified: Feb. 17, 2017
%-----------------------------------------------------------------------
 
tcSENSE  = tic;

A        = param.A;
At       = param.At;
U        = param.U;
Ut       = param.Ut;

iIter    = param.iIter; % Inner iteration
oIter    = param.oIter; % Total Outer iteration
minIter  = param.minIter;
stpThrsh = param.stpThrsh;

L1       = param.L1;

lmb0     = param.lmb0;
sFtri    = param.sFtri;
sFtrf    = param.sFtrf;
fstIter  = param.fstIter;

lRes     = param.lRes; 

%Use the power iteration to estimate the Lipschitz constant if not provided
if isempty(L1)
    L1 = 2.05*powerIter(param.A, param.At, size(x0));
end

% Let sFtr start with a high value; this help with high snr (snr>60 datasets);
sFtr = [sFtri*10.^(log10(sFtrf/sFtri) .* ((1:round(oIter*0.75))'/round(oIter*0.75)).^4); sFtrf*ones(oIter - round(0.75*oIter),1)];

x = x0; clear x0;
z = x;
t = 1;


[objDF, objWav] = wObjective(b, x, lmb0, param);
disp(['Initial objective function: ' sprintf('%0.3f',objDF) ' + ' ...
      sprintf('%0.3f',objWav) ' = ' sprintf('%0.3f',objDF+objWav)]);

dispim = figure;
for o = 1:oIter % outer (reweighting) iteration
    stop = 0;
    iter = 0;
    
    % caculate regularization weights for each subband
    Ux = U(x);
    param.sFtr = sFtr(o);
    param.Nlim = lRes(2);
    if o>=ceil(lRes(1)*oIter)
        param.Nlim = inf;
    end
    [lmb, ~] = findWFst(Ux, param);
    clear Ux;
    
    % For the first outer iteration
    if oIter == 1
        fI = fstIter;
    else
        fI = 1;
    end
    
    
    while iter<=(fI*iIter) && stop==0 % inner iteration that solves the convex problem
        iter = iter+1;
        
        coefG = U(z - 1/L1 * wGradient(b, z, A, At));
        coefG = shrink(coefG, lmb/L1);
        x_new = Ut(coefG);
        t_new = (1 + sqrt(1 + 4*t^2))/2;
        z = x_new + (t - 1)/t_new*(x_new - x);
        x_old = x;
        x = x_new;
        t = t_new;
        
        if iter>=(fI*minIter*iIter)
            tmp = (norm(x_new(:) - x_old(:)) / norm(x_new(:)));
            if tmp  < stpThrsh
                disp('Stopping threshold reached, terminating inner loop');
                stop = 1;
            end
        end
        
    end
    
    [objDF, objWav] = wObjective(b, x, lmb, param);
    tmp = (norm(x_new(:) - x_old(:)) / norm(x_new(:)));
    
    %2D cine image (show the image of central-frame)
    im0 = x(:,:,ceil(size(x,3)/2),1,1,1);
    figure(dispim);
    imagesc(abs(im0),[0,0.3*max(abs(x(:)))]); axis('image'); colormap(gray); title(['oIter: ' num2str(o)]);
    pause(0.01);
    disp(['OuterIter: ' num2str(o)  ', Change in x: ' sprintf('%0.3e',tmp) ', objective: ' sprintf('%0.3f',objDF) ' + ' sprintf('%0.3f',objWav) ' = ' sprintf('%0.3f',objDF+objWav)]);
    
end
disp(['Total recovery time: ' sprintf('%0.2f', toc(tcSENSE)) ]);
