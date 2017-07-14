function omega1_rms = compute_omega1_rms(pulse)
%{
omega1_rms = compute_omega1_rms(pulse) function computes mean omega1 value 
of shaped MT pulse using Yarnykh and Yuan model: MT shaped pulse of 
duration t_mt  approximated by cw rectangular pulse of the same duration 
and equivalent average power

input:
    pulse - data structure containig all MT pulse properties
output:
    omega1_rms - averaged omega1 value, [rad/s]
%}

t_mt = pulse.t_mt;
omega1_sq_func = pulse.omega1_sq_func;  % function handle describing power 
                                        % of MT pulse depending on time
integral_value = integral(omega1_sq_func, 0, t_mt);
omega1_rms = sqrt(integral_value / t_mt) * pulse.amp * 2 *pi;

end

