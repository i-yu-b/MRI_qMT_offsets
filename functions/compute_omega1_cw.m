function omega1_cw = compute_omega1_cw(pulse)
%{
omega1_cw = compute_omega1_cw(pulse) function computes mean omega1 value 
of shaped  MT pulse using Sled and Pike CW approach: MT shaped pulse of 
duration t_mt approximated by cw rectangular pulse of TR duration and 
equivalent average power as shaped MT pulse
input:
    pulse - data structure containig all MT pulse properties
output:
    omega1_cw - averaged omega1 value, [rad/s]
%}

TR = pulse.TR;
t_mt = pulse.t_mt;
omega1_sq_func = pulse.omega1_sq_func;  % function handle describing power 
                                        % of MT pulse depending on time
integral_value = integral(omega1_sq_func, 0, t_mt);
omega1_cw = sqrt(integral_value / TR) * pulse.amp * 2 * pi;

end
