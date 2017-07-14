function Mzf = mtoff_fit(parameters,offsets,pulse,omega1_method, ...
                        Mz_method,lineshape,sig)
%{
function Mz = mtoff_fit(parameters,constants) solves 2-compartment
bloch equations in the presence of off-resonance saturation pulses 

input:
    parameters:
       parameters(1) - macromolecular-to-free pool exchange rate,[s-1]
       parameters(2) - R1 relaxation rate of macromolecular pool, 
                        commonly turned off as a fitting parameter and 
                        assumed to be 1 s-1
       parameters(3) - R2 relaxation rate of macromolecular pool, [s-1]
       parameters(4) = kmf*Mom/R1f,where
                                  Mom - Mo value of macromolecular pool, 
                                  R1f - R1 relaxation rate of free pool
       parameters(5) = 1/(R1f*T2f), where 
                                     T2f - T2 relaxation time of free pool
    
    offsets - array of frequency offsets for MT pulse, [Hz]
    pulse - data structure containig all MT pulse properties
    
    omega1_method - function handle for function calculating mean omega1
                    values for different MT pulse approximations
                    currently available: 
                    @compute_omega1_cw - Sled and Pike continious wave                                         
                        (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)

                    @compute_omega1_rp - Sled and Pike rectangular pulse                                        
                        (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)

                   @compute_omega1_rms - Yarnykh and Yuan rectangular pulse 
                          (Yarnykh VL, Yuan C. 2004. Neuroimage 23:409-424)
    
    Mz_method - function handle for function calculating Mz value for 
                different MT pulse approximations
                currently available: 
                  @SPGR_cw - Sled and Pike continious wave model 
                  (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)

                  @SPGR_rp - Sled and Pike rectangular pulse model
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
%}

% calculate integral under macromolecular pool spectrum
T2m = parameters(3); %T2 relaxation time of macromolecular pool

G = lineshape(offsets, T2m);

% calculate averaged omega1 value
omega1_aver = omega1_method(pulse);
% calculate mean saturation rate of macromolecular pool
W = G .* pi .* omega1_aver.^2;

Mzf =  Mz_method(parameters, W, offsets, omega1_aver);
semilogx(offsets, sig, '-o',offsets, Mzf,'-k');
drawnow

