function mtoffSet = mtoffAnalysis(dset)
%{ 
function mtoffSet = mtoffAnalysis(dset) performs MT offset saturation 
analysis on the dataset

inputs:
           dset.img - 2D or 3D image data in the format:
                     (x, y, z(optional), offsets, flip_angles)
           dset.mask - optional mask for processing data 
           dset.pars - set of all necessary parameters fo data analysis
outputs:  mtoffSet - a data set containing main MT offset saturation 
                      parameter maps of: 
           mtoffSet.kfm - exchange rate free-to-macromolecular pool,[s-1]
           mtoffSet.kmf - exchange rate macromolecular-to-free pool,[s-1]
           mtoffSet.R1f - relaxation rate of free pool, [s-1]
           mtoffSet.PSR - pool-size ratio
           mtoffSet.T2m - relaxation time of macromolecular pool, [s-1]
           mtoffSet.T2f - relaxation time of free pool, [s-1]
%}

% load in the dataset
sz = size(dset.img); 

% load in the needed parameters
gamma = 42.577478; % MHz/T

tr = dset.pars.methpars.PVM_RepetitionTime;
mtoffSet.kfm = zeros(size(mask));s

offsets = dset.pars.methpars.PVM_MagTransFL; % offset freq in [Hz]
flip_angles = dset.pars.methpars.PVM_MagTransPower; % B1 values in [uT]
pulse_shape = dset.pars.methpars.PVM_MagTransPulse1Enum; 
% other pulse parameters
% load B0 map
% load B1 map

% load T1 map in [s]
T1_map = t1.T1/1000;

% choose lineshape
% choose method

% build pulse

% define a mask if one is not given
if isfield(dset,'mask')
    mask = dset.mask;
else
    mask = true(prod(sz(1:3),1));
end

% initialize the mtoff dataset
mtoffSet.kfm = zeros(size(mask));
mtoffSet.kmf = zeros(size(mask));
mtoffSet.PSR = zeros(size(mask));
mtoffSet.BPF = zeros(size(mask));
mtoffSet.R1f = zeros(size(mask));
mtoffSet.T2m = zeros(size(mask));
mtoffSet.T2f = zeros(size(mask));

lb = [0 0 0 0 0];
ub = [inf inf inf inf inf];

tot_evals = sum(mask(:));
evals = 0;

warning('off','MATLAB:singularMatrix')

fprintf('%3.0f %% done...',0);
for ro=1:size(dset.img,1)
    for pe=1:size(dset.img,2)
        for sl=1:size(dset.img,3)
            if mask(ro,pe,sl)
                
                sig = squeeze(abs(dset.img(ro,pe,sl,:)));
                % normalize signal
                sig = sig/max(sig);

                % initial guess
                b0 = [70,5,8e-6,10,20]; 
                
                % fit the data
                opts = optimset('display','off');
                [b,~,res,~,~,~,jac] = lsqnonlin(@(x) ...
                         mtoff_fit(x,"all parameters")-sig,b0,lb,ub,opts);
     
                % extract map parameters
                mtoffSet.kmf(ro,pe,sl) = b(1);
                % relaxation rate of macromolecular pool has little
                % influence on fitting and typically assumed to be 1s-1
                % alternatively can fit out of experimental data = b(2)
                R1m = 1; 
                
                R1_obs = 1/T1_map(ro,pe,sl);
                mtoffSet.R1f(ro,pe,sl) = R1_obs/...
                               (1 + b(4)*(R1m - R1_obs)/(R1m-R1_obs+b(1)));
                
                mtoffSet.kfm(ro,pe,sl) = b(4)*mtoffSet.R1f(ro,pe,sl);                   
                mtoffSet.PSR(ro,pe,sl) = mtoffSet.kfm(ro,pe,sl)/b(1);
                mtoffSet.BPF(ro,pe,sl) = mtoffSet.PSR(ro,pe,sl)/...
                                         (1+ mtoffSet.PSR(ro,pe,sl));
                mtoffSet.T2m(ro,pe,sl) = b(3);
                mtoffSet.T2f(ro,pe,sl) = 1/(b(5)*R1f);

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


