function Mzf = SPGR_ramani(parameters, W, offsets, omega1)
%{
function Mz = SPGR_ramani(parameters, W) calculates Mz value using Ramani's
model (Ramani A, et al. 2002. Magn Reson Imaging 20:721-731).

input:
    parameters:
       parameters(1) - kmf, macromolecular-to-free pool exchange rate,[s-1]
       parameters(2) - R1 relaxation rate of macromolecular pool, 
                        commonly turned of as a fitting parameter and 
                        assumed to be 1 s-1
       parameters(3) - R2 relaxation rate of macromolecular pool, [s-1]
       parameters(4) = R*Mom/R1f,where
                                  Mom - Mo value of macromolecular pool, 
                                  R1f - R1 relaxation rate of free pool
       parameters(5) = 1/(R1f*T2f), where 
                                     T2f - T2 relaxation time of free pool
    W - averaged saturation rate of macromolecular pool, [s-1]
    omega1 - approximated omega1 value using Ramani's model
    offsets - array of frequency offsets for MT pulse, [Hz]

output:
    Mzf - normalized(!) proton magnetization vector
%}

kmf = parameters(1); 
R1m = 1; % common assumption. To fit out the parameter uncomment 
         % the following line:
% R1m = parameters(2);
T2m = parameters(3); 
coeff1 = parameters(4); % kmf*Mom/R1f
coeff2 = parameters(5); % 1/(R1f*T2f)

% calculating saturation rate of free pool * R1f
Wf = (omega1 ./ (2 * pi * offsets)).^2 * coeff2; 

numerator = R1m * coeff1 + W + R1m + kmf; 
denominator = coeff1 * (R1m + W) + (1 + Wf).*(R1m + W + R);
Mzf = numerator./denominator;

end

