function Mzf_norm = SPGR_yarnykh(parameters, Wm, offsets, omega1, R1f_obs, pulse)
%{
function Mz = SPGR_yarnykh(parameters, W) calculates Mz value using Yarnykh
              solution for MT off-res saturation
              (Yarnykh VL, Yuan C. 2004. Neuroimage).
input:
    parameters:
    parameters(1) - pool-size ratio
    parameters(2) - kmf - exchange rate macromolecular-to-free pool,[s-1]
    parameters(3) -  relaxation rate of free pool, [s-1]
    parameters(4) - relaxation rate of macromolecular pool, [s-1]
    parameters(5) - T2 relaxation time of free pool, [s-1]
    parameters(6) - T2 relaxation time of macromolecular pool, [s-1]
                                     
    Wm - averaged saturation rate of macromolecular pool, [s-1]
    omega1 - array of averaged omega1 values, []
    offsets - array of frequency offsets for MT pulse, [Hz]
    R1f_obs - measured from separate experiment R1 relaxation of free pool,
              [s-1]

output:
    Mzf - normalized(!) proton magnetization vector
%}

PSR   = parameters(1);
kmf  = parameters(2);
% R1f = parameters(3); uncomment if R1 map is not provided
R1m = parameters(4);
T2f = parameters(5);
T2m = parameters(6);
kfm = kmf * PSR;
BPF   =  PSR / (1+PSR);       % Bound Pool Fraction

% common assumption:
% R1m = 1; % or:
% R1m = R1f; % uncomment if needed
%{ 
the latter is preferred in Yarnykh, because "it eliminates the
question about magnetic field dependence of RB1 and can
be uniformly applied at different field strengths."
%}

R1f = R1f_obs - kfm *(R1m - R1f_obs) / (R1m - R1f_obs + kfm/PSR);

% extract additional time parameters from pulse datastructure
% Kevin, add this structure fields in function build_pulse

TR = pulse.TR;
tm  =  pulse.t_mt;       % saturation offset pulse duration
ts  =  pulse.ts;       % delay before excitation pulse
tr  =  pulse.tr;       % delay after the excitation pulse
alpha_ex = pulse.angle /180 * pi;   % excitation pulse flip angle, [rad]


% the vector of equilibrium magnetization with elements Meq = [1, PSR]';
% Attention, it's different from original paper: here data is normalized to
% Meq value of free water, in the paper it's normalized to total magnetiza-
% tion. Meaning also, that possibly Cabana code is incorrect.

Meq = [1, PSR]';
                     
% rotation of the magnetization by the excitation pulse with the flip angle 
% alpha_ex
C = diag([cos(alpha_ex) 1]);
% unit matrix
I   =  eye(2);
% longitudal relaxation matrix
Rl  =  [ -R1f-kfm,  kmf;  kfm,  -R1m-kmf ];
% longitudinal relaxation during delays before and after excitation pulse
Es  =  expm(Rl*ts);
Er  =  expm(Rl*tr);
% averaged absorption rate of free pool, [s-1]
Wf = (omega1 ./ (2 * pi * offsets)).^2 * (1 / T2f); 

%{
 Kevin, in case of Ramani approach we used this equation:
Wf = (omega1 ./ (2 * pi * offsets)).^2 * (1 / T2f); 
In Cabana, 2016 paper they suggest using this formula for Yarnykh method 
as well. But it's correct only in case of relatively large offsets >> spectral
width of free pool. Stricter way is to calculate it properly:
WF = pi * omega1_averaged * G(offsets, T2f), G - superlorentzian lineshape.
%}
% convinient value for calculing Mz steady state vector
A   =  R1f*R1m + R1f*kmf + R1m*kfm;

% loop over different powers of saturation pulse and offsets. Maybe
% flattening W to an array vs keeping it as a matrix to get rid of one inner
% loop?
num_powers = length(omega1);
num_offsets = length(offsets);

Mzf = zeros(num_offsets, num_powers);
for i=1:num_powers
    for j=1:num_offsets
        % saturation matrix 
        W   =  -diag([Wf(j, i) Wm(j, i)]);       
        % convinient value for calculing Mz steady state vector
        D   =  A + (R1f+kfm)*Wm(j,i) + (R1m+kmf)*Wf(j,i) + Wm(j,i)*Wf(j,i); 
    % Mz steady state vector 
    % if you want to choose normalization to 
    % total magnetization need to use [(1-BPF); BPF] coefficients instead
    % of [1, PSR].
        Mzss =  1./D *[A + R1f*Wm(j, i); ...
                   PSR*(A + R1m*Wf(j, i))];   
     
        Em  =  expm( (Rl+W)*tm );   
        denominator  =  I - Es * Em * Er * C;
        numerator  =  (Es*Em*(I-Er) + (I-Es))*Meq + Es*(I-Em)*Mzss;
        Mz =  numerator \ denominator;
        Mzf(j,i) = Mz(1);
    end
end

% Normalize to the same sequence without MT prepulse to exclude factors 
% associated with M0,T2*, coil sensitivity.
% Wf = 0; Wm = 0

    Em0 = expm(Rl*tm);
    D0 = A;
    Meq0 =  [ 1; PSR];
    Mzss0 =  1/D0 * [A; PSR*A];
    denominator0  =  I - Es * Em0 * Er * C;
    numerator0  =  (Es*Em0*(I-Er) + (I-Es))*Meq0 + Es*(I-Em0)*Mzss0;
    Mz0 =  numerator0 \ denominator0;
    Mzf0 = Mz0(1);
   
    Mzf_norm = Mzf ./ Mzf0;
end
