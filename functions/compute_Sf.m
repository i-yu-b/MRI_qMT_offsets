function Sf = compute_Sf(offsets, T2f, pulse)
%{
function Sf = compute_Sf(flip_angles,offsets,T2f)
         computes Sf value describing partial saturation effect of MT pulse
         on free pool. For details see:
         (Sled JG, Pike GB. 2001. Magn Reson Med 46:923-931) and
         (Pike GB. 1996. Magn Reson Med 36:95?103.)
input:
        offsets - array of frequency offsets for MT pulse, [Hz]
        flip_angles - array of flip angles of MT pulse, [rad] -
        KEVIN, some people dont read out amp of pulses from Paravision,
        but rather calculate if from flip angle, don't know which one is 
        better.

        T2f - relaxation time of free pool, [s-1]
        pulse - data structure describing properties of MT pulse
output:
        Sf - saturation fraction of free water pool
%}

amps = pulse.amp;
pulse_shape = pulse.omega1_func;
M0f = [0 0 1];
Sf = zeros(length(offsets),length(amps));

for i=1:length(offsets)
    for j=1:length(amps)
        [~, Mf_sol] = ode23(@(t,Mf) ...
                Bloch_eq(t, Mf, T2f, pulse_shape, offsets(i), amps(j)), ...
                [0 pulse.t_mt], M0f);
        Sf(i,j) = Mf_sol(end,3);
    end
end
end
