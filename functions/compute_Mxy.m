function Mxyf = compute_Mxy(kmf, R1m, kfm, R1f, Sf, Wr, TR, flip_angles)
%{
function Mxyf = compute_Mxy(kmf, R1m, kfm, R1f, Sf, Wr, TR, flip_angles)
         computes Mxy value using signal equation based on Sled and Pike 
         continious wave model 
         (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)
         Kevin, you also can include Cabana paper (qMT lab) 2016, as they
         have nicely written equations.
input:
    parameters:
    kmf - exchange rate macromolecular-to-free pool,[s-1]
    R1m - relaxation rate of macromolecular pool, [s-1]
    kfm - exchange rate free-to-macromolecular pool,[s-1]
    R1f - relaxation rate of free pool, [s-1]
    Sf - correction for saturation effect of MT pulse on free pool
    Wr - averaged saturation rate of macromolecular pool, [s-1]
    TR - repetition time of pulse sequence, [s]
    flip_angles - array of flip angles of MT pulse, [rad]

output:
    Mxyf - Mxy proton magnetization vector
%}

% calculate Mzf steady state (assuming Mz0 = 1)
Mz_numerator =  R1m * kfm + R1m * R1f + R1f*kmf + Wr*R1f;
Mz_denominator  =  R1m * R1f + R1m * kfm + R1f * kmf + Wr * R1f + Wr * kfm;
Mzf   =  Mz_numerator ./ Mz_denominator;

sqrt_value  =  sqrt( (R1f + kfm + R1m + kmf + Wr).^2 - ...
                      4*(R1f* R1m + kfm*R1m + kmf*R1f + R1f*Wr + kfm*Wr));
lambda1  =  (R1f + kfm + R1m + kmf + Wr + sqrt_value)/2;
lambda2  =  (R1f + kfm + R1m + kmf + Wr - sqrt_value)/2;

E1  =  exp( -lambda1*TR );
E2  =  exp( -lambda2*TR );

numerator  =  (E1-1) .* (E2-1) .* (lambda2-lambda1) .* ...
               Sf .* Mzf .* sin(flip_angles);
           
denominator  =  (E1-1) .* (Sf.*E2-1) .* (lambda2-lambda1) + ...
                (Sf-1) .* (E2-E1) .* (lambda2 - R1f - kfm);
Mxyf =  numerator ./ denominator;
end

