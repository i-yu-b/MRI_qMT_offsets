function Mxyf = compute_Mxy_rp(kmf, R1m, kfm, R1f, S, Wm, TR, flip_angle, tau, M0)
%{
function Mxyf = compute_Mxy_rp(kmf, R1m, kfm, R1f, Sf, Wr, TR, flip_angles)
         computes Mxy value using signal equation based on Sled and Pike 
         rectangular pulse model 
         (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931)
         
input:
    parameters:
    kmf - exchange rate macromolecular-to-free pool,[s-1]
    R1m - relaxation rate of macromolecular pool, [s-1]
    kfm - exchange rate free-to-macromolecular pool,[s-1]
    R1f - relaxation rate of free pool, [s-1]
    Sf - correction for saturation effect of MT pulse on free pool
    Wm - averaged saturation rate of macromolecular pool, [s-1]
    TR - repetition time of pulse sequence, [s]
    flip_angle - flip angle of MT pulse, [rad]
    M0 = [1; PSR] 
    tau - full-width at half-maximum value of MT shaped pulse, [s-1]

output:
    Mxyf - Mxy proton magnetization value
%}

    M0f = M0(1);
    M0m = M0(2);
    
    % free precession
    Afp  = [R1f+kfm,  -kmf; -kfm,  R1m+kmf];
    expAfp = expm(-(TR-tau)*Afp);
    
    % continious saturation
    Acw  =  [R1f+kfm,  -kmf ; -kfm,  R1m+kmf+Wm];
    expAcw = expm( -tau/2.*Acw);
    
    % steady state of the magnetization established after a long period of 
    % continuous-wave irradiation of the macromolecular pool
    denominator =  R1m  *R1f + R1m * kfm + R1f * kmf + Wm * R1f + Wm * kfm;
    numerator_f = M0f * (R1m * kfm + R1m * R1f + R1f * kmf + Wm * R1f);
    numerator_m = M0m *(R1m * R1f + R1m * kfm + R1f * kmf);
    Mzf_ss = numerator_f / denominator; 
    Mzm_ss = numerator_m/ denominator; 
    
    Mz_ss = [Mzf_ss; Mzm_ss];
    I = eye(2);

%     Mz = (I - expAcw * expAfp * expAcw * diag([Sf.*cos(alpha) Sr]))\ ...
%         ( (I + expAcw * (-I + expAfp * (I -  expAcw)))*Mss + expAcw*(I - expAfp)*M0_inf);
    
    Mxyf = Mz(1,:).*sin(flip_angle).*S(:,1);

end