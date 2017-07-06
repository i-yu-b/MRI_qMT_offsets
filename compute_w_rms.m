function w_rms = compute_w_rms(pulse)
%{
w_rms = compute_w_rms(pulse) function computes mean omega1 value of shaped 
MT pulse using Yarnykh and Yuan model: MT shaped pulse of duration t_mt 
approximated by cw rectangular pulse of the same duration and equivalent 
average power

input:
    pulse - data structure containig all MT pulse properties
output:
    w_rms - averaged omega1 value
%}

t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);
w_rms = sqrt(integral_value / t_mt);

end

