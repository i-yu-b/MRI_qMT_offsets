function Mxyf = SPGR_rp(parameters, Wm, offsets, omega1, R1f_obs, tau)
%{
function Mxyf = SPGR_rp(parameters, Wr, offsets, omega1, R1f_obs) 
         calculates Mxy value using Sled and Pike rectangular pulse model 
         solution for MT off-res saturation
         (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)
              
input:
    parameters:
    parameters(1) - pool-size ratio
    parameters(2) - kmf - exchange rate macromolecular-to-free pool,[s-1]
    parameters(3) -  relaxation rate of free pool, [s-1]
    parameters(4) - relaxation rate of macromolecular pool, [s-1]
    parameters(5) - relaxation time of free pool, [s-1]
    parameters(6) - relaxation time of macromolecular pool, [s-1]                                  
    Wm - averaged saturation rate of macromolecular pool, [s-1]
    omega1 - array of averaged omega1 value, []
    offsets - array of frequency offsets for MT pulse, [Hz]
    R1f_obs - measured from separate experiment R1 relaxation of free pool,
              [s-1]
    tau - full-width at half-maximum value of MT shaped pulse, [s-1]
output:
    Mxyf - normalized(!) proton magnetization vector
%}

    PSR   = parameters(1);
    kmf  = parameters(2);
% R1f = parameters(3); uncomment if R1 map is not provided
    R1m = parameters(4);
    T2f = parameters(5);
    T2m = parameters(6);
    kfm = kmf * PSR;

% common assumption:
% R1m = R1f; % uncomment if needed

    R1f = R1f_obs - kfm *(R1m - R1f_obs) / (R1m - R1f_obs + kfm/PSR);

    TR = pulse.TR;
    flip_angles = pulse.flip_angle /180 * pi;   % MT pulse flip angles, [rad]

% saturation of the free pool due to the MT pulse. 
% Kevin, you  speed up fitting by make Sf pre-computed, because it only
% depends on shape of pulse, its amplitude, power and offset. So you can 
% create table Sf(offset, power).
    Sf = compute_Sf(flip_angles,offsets,T2f);
    
 %{
    Kevin, according to rp approximation pulse sequence has the following 
    steps: instantaneous saturation of the free pool from the MT and 
    excitation pulse, continuous-wave irradiation of the restricted
    pool for a period tau/2, a period TR - tau of free precession,
    and finally another period of continuous-wave irradiation of
    duration tau/2.
 
    Therefore you need to pass tau value - full-width at half-maximum value
    of MT shaped pulse. It's better to calculate in the beginning, 
    because we also use this value to calculate averaged power. You 
    can use fwhm.m function I've sent you
   
  %}
    
    % now have to iterate across different powers (compared to Sled&Pike)
    num_powers = length(omega1);
    M0 = [1; PSR];
    Sr = 1; % alternatively can be changed
    S = [Sf; Sr];
    
    % loop over different powers of saturation pulse and offsets. Maybe
    % flattening W to an array vs keeping it as a matrix to get rid of one inner
    % loop?
    num_powers = length(omega1);
    num_offsets = length(offsets);

    Mxy = zeros(num_offsets, num_powers);
    for i=1:num_powers
        for j=1:num_offsets
            Mxy(j, i)  = compute_Mxy_rp(kmf, R1m, kfm, R1f, S,  Wm(j,i), ...
                                        TR,flip_angles(i), tau, M0);                              
        end
    end
    % no MT pulse (for normalization) Sf = 1; Wm = 0; Sr = 1;
    Mxy0  = compute_Mxy_rp(kmf, R1m, kfm, R1f, [1; 1], 0, TR, ...
                           flip_angles, tau, M0);
                             
    Mxyf = Mxy ./ Mxy0; % normalized value

end


