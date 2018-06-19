%%
%Input:
%kdata: data in k-space
%       [kx, ky, kz, coil, time] or [kx,ky,coil,time] (other dims such as "cha","slc" must be 1)
%samp: sampling pattern, same for all coils, so no coil information
%      it is fully sampled in RO direction always.
%      two dimensions less than kdata
%      [kx, ky, kz, time] or [kx ky time]
%      type: boolean, true indicating the point was acquired
%param:
%     .mthd  : method
%     .fil   : order of smoothing matrix, default 9
%     .opt   : type of phase correction

%Output:
%S: sensitivity map
%cI: combined Image
% Chong Chen @OSU 06/18/2018

%%
function [S, cI] = coilSen(kdata, samp, param)
% works both for 2D and 3D

kdata = squeeze(kdata);
samp = squeeze(samp);
%% 2D to pesudo 3D
if ndims(kdata) == 4 %kx,ky,coil,frame
    kdata = permute(kdata,[1,2,5,3,4]); %kx,ky,kz,coil,frame
    samp = permute(samp,[1,2,4,3]); %kx,ky,kz,frame
end

if ndims(kdata) == 5 % kx, ky, kz, coil, time
    % average over time for each coil
    kdata_avg = sum(kdata,5)./(repmat(sum(samp,4),[1,1,1,size(kdata,4)])+eps);
    samp_avg = logical(sum(samp,4));
elseif ndims(kdata) == 4
    kdata_avg = kdata;
end
    disp(['Espirit,' num2str(1) ' sensitivity maps']);
    CalibSize = size(samp_avg);
    if numel(CalibSize) == 2, CalibSize(3) = 1;end
    CalibSize = min(CalibSize,param.ACSsz); % specify the calibration region
    [cI,S] = espirit_sens3D(kdata_avg,samp_avg,CalibSize,param);
    [cI,S] = sensCorrect3D(cI, S, param.avgPhs); % Remove phase from time average    
    S = squeeze(S); cI = squeeze(cI); % pesudo-3D to 2D






