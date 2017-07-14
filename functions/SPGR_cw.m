function Mxyf = SPGR_cw(parameters, Wr, offsets, omega1, R1f_obs)
%{
function Mxyf = SPGR_cw(parameters, Wr, offsets, omega1, R1f_obs) 
         calculates Mxy value using Sled and Pike continious wave model 
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
    Wr - averaged saturation rate of macromolecular pool, [s-1]
    omega1 - averaged omega1 value, []
    offsets - array of frequency offsets for MT pulse, [Hz]
    R1f_obs - measured from separate experiment R1 relaxation of free pool,
              [s-1]

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

% extract additional time parameters from pulse datastructure
% Kevin, add {alpha_ex} structure fields in function build_pulse

TR = pulse.TR;
flip_angles = pulse.flip_angle /180 * pi;   % MT pulse flip angles, [rad]

% saturation of the free pool due to the MT pulse. 
% Kevin, you  speed up fitting by make Sf pre-computed, because it only
% depends on shape of pulse, its amplitude, power and offset. So you can 
% create table Sf(offset, power).
Sf = compute_Sf(flip_angles,offsets,T2f);

Mxy  = compute_Mxy(kmf, R1m, kfm, R1f, Sf, Wr, TR, flip_angles);
% no MT pulse (for normalization)
Mxy0 = compute_Mxy(kmf, R1m, kfm, R1f, 1, 0, TR, flip_angles); 

Mxyf = Mxy ./ Mxy0; % normalized value

end