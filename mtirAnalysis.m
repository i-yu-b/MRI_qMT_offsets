function mtoffSet = mtoffAnalysis(dset)
% function mtoffSet = mtoffAnalysis(dset) performs MT offset saturation 
% analysis on the dataset
%
% inputs:
%           dset.img - 2D or 3D image data in the format:
%                     (x, y, z(optional), offsets, flip_angles)
%           dset.mask - optional mask for processing data 
%           dset.pars - set of all necessary parameters fo data analysis
% outputs:  mtoffSet - a data set containing main MT offset saturation 
%                      parameter maps of: 
%           mtoffSet.R - exchange rate, s-1
%           mtoffSet.kmf - exchange rate macromolecular-to-free pool, s-1
%           mtoffSet.R1f - relaxation rate of free pool, s-1
%           mtoffSet.PSR - pool-size ratio
%           mtoffSet.T2m - relaxation time of macromolecular pool, s-1
%           mtoffSet.T2f - relaxation time of free pool, s-1


% load in the dataset
sz = size(dset.img); 

% load in the needed parameters
tr = dset.pars.methpars.PVM_RepetitionTime

offsets = dset.pars.methpars.PVM_MagTransFL; % offset freq in Hz
flip_angles = dset.pars.methpars.PVM_MagTransPower; % B1 values in uT
pulse_shape = dset.pars.methpars.PVM_MagTransPulse1Enum 
% other pulse parameters
% load B0 map
% load B1 map
gamma = 42.577478 % MHz/T

% define a mask if one is not given
if isfield(dset,'mask')
    mask = dset.mask;
else
    mask = true(prod(sz(1:3),1));
end

% initialize the mtoff dataset
mtoffSet.M0a = zeros(size(mask));
mtoffSet.M0b = zeros(size(mask));
mtoffSet.PSR = zeros(size(mask));
mtoffSet.BPF = zeros(size(mask));
mtoffSet.kmf = zeros(size(mask));
mtoffSet.T1 = zeros(size(mask));
mtoffSet.ci = cell(size(mask));

lb = [0 0 0 0 -1];
ub = [inf inf inf inf 1];

tot_evals = sum(mask(:));
evals = 0;

warning('off','MATLAB:singularMatrix')

fprintf('%3.0f %% done...',0);
for ro=1:size(dset.img,1)
    for pe=1:size(dset.img,2)
        for sl=1:size(dset.img,3)
            if mask(ro,pe,sl)
                
                sig = squeeze(abs(dset.img(ro,pe,sl,:)));

                % initial guess for soft IR, alpha = 0.8
                b0 = [max(sig)/2, max(sig)/3,10,2,-0.9]; 
                
                % fit the data
                opts = optimset('display','off');
                [b,~,res,~,~,~,jac] = lsqnonlin(@(x) remmi.util.sir(x,ti',td)-sig,b0,lb,ub,opts);
                
                % load the dataset
                mtirSet.M0a(ro,pe,sl)=b(1);
                mtirSet.M0b(ro,pe,sl)=b(2);
                mtirSet.PSR(ro,pe,sl)=b(2)/b(1); %M0b/M0a
                mtirSet.BPF(ro,pe,sl)=b(2)/(b(2)+b(1));
                mtirSet.kmf(ro,pe,sl)=b(3);
                mtirSet.T1(ro,pe,sl)=1/b(4);
                
                % save confidence intervals on the original parameters
                mtirSet.ci{ro,pe,sl} = nlparci(b,res,'jacobian',jac); 
                
                evals = evals+1;
            end
        end
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f %% done...',evals/tot_evals*100);
    end
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%3.0f %% done...\n',100);

warning('on','MATLAB:singularMatrix')


