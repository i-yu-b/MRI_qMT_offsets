function mtoffSet = mtoffAnalysis(dset, method,lineshape)
%{ 
function mtoffSet = mtoffAnalysis(dset) performs MT offset saturation 
analysis on the dataset

inputs:
           dset.img - 2D or 3D image data in the format:
                     (x, y, z(optional), offsets, flip_angles)
           dset.mask - optional mask for processing data 
           dset.pars - set of all necessary parameters fo data analysis
           
           method_name - chosen method for performing off res MT analysis 
                         shaped pulse. Choose between:
                         'sled_cw' - Sled and Pike continious wave model                                         
                                   (Sled JG, Pike GB. 2001. Magn Reson Med)
                         'sled_rp' - Sled and Pike rectangular pulse model                                       
                                (   Sled JG, Pike GB. 2001. Magn Reson Med)
                         'ramani' - Ramani model
                                (Ramani A, et al. 2002. Magn Reson Imaging)
                                    (Yarnykh VL, Yuan C. 2004. Neuroimage)   
                         'yarnykh' - Yarnykh and Yuan rectangular pulse 
                                (Yarnykh VL, Yuan C. 2004. Neuroimage)
                         By default it's 'ramani'
           Depending on method_name, two following parameters are set up:

           omega1_method_name - chosen method for approximation of MT 
                                shaped pulse. Currently available:  
                       = @compute_omega1_cw - Sled and Pike continious wave                                         
                       = @compute_omega1_rp - Sled and Pike rectangular pulse                                        
                       = @compute_omega1_rms - Yarnykh and Yuan rectangular pulse 
          

            Mz_method - function handle for function calculating Mz value 
                        for different MT pulse approximations.
                        Currently available: 
                  @SPGR_cw - Sled and Pike continious wave solution 
                  (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)

                  @SPGR_rp - Sled and Pike rectangular pulse solution
                  (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)

                  @SPGR_ramani - Ramani solution for MT off-res saturation
                  (Ramani A, et al. 2002. Magn Reson Imaging 20:721-731)
                  
                  @SPGR_yarnykh - Yarnykh solution for MT off-res saturation
                  (Yarnykh VL, Yuan C. 2004. Neuroimage)
                        
    lineshape - function handle for function calculating lineshape of
                macromolecular pool. 
                currently available:
                  @gaussian - for gaussian lineshape
                  @lorentzian - for lorentzian lineshape
                  @superlorentzian - for superlorentzian lineshape
                  @superlorentzian_res - for superlorentzian lineshape with
                                         spline of values for offsets less
                                         than 1.5 kHz.
outputs:   mtoffSet - a data set containing main MT offset saturation 
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

    % set up default method
    if ~exist('method','var')
     % set default
     method = 'ramani';
    end
    
    % set up methods for calculating Mz and averaged omega1:
    if strcmp(method, 'sled_cw')
        Mz_method = @SPGR_cw;
        omega1_method = @compute_omega1_cw;
    elseif strcmp(method, 'sled_rp')
        Mz_method = @SPGR_rp;
        omega1_method = @compute_omega1_rp;
    elseif strcmp(method, 'ramani')
        Mz_method = @SPGR_ramani;
        omega1_method = @compute_omega1_cw;
    elseif strcmp(method, 'yarnykh')
        Mz_method = @SPGR_yarnykh;
        omega1_method = @compute_omega1_rms;
        
    % load in the needed parameters
    gamma = 42.577478; % MHz/T
    offsets = dset.pars.methpars(1).PVM_MagTransFL'; % offset freq in [Hz]
                                                     % column vector
    
    % build pulse data structure
    pulse = build_pulse(dset.pars);
    
    % define a mask if one is not given
    if isfield(dset,'mask')
        mask = dset.mask;
    else
        mask = true(prod(sz(1:3),1));
    end
    
    % load B0 map

    % load B1 map

    % load T1 map in [s]
    T1_map = t1.T1/1000;

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
                
                sig = squeeze(abs(dset.img(ro,pe,sl,:,:)));
                % normalize signal for every power
                sig = sig./max(sig,[],1);

                % initial guess
                b0 = [70,5,8e-6,10,2]; 
                lb = [0, 0, 0, 0, 0];
                ub = [inf, inf, inf, inf, inf];
            
                
                % fit the data
                opts = optimset('display','off');
  %{
    needs to be re-written as different methods have different parameters
    to fit! Currently set up for Ramani method.            
    %}
               [b,~,res,~,~,~,jac] = lsqnonlin(@(x) ...
                     mtoff_fit(x,offsets, pulse, @compute_omega1_rms,...
                     @SPGR_ramani,@superlorentzian,sig)-sig,b0,lb,ub,opts);

                    
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


