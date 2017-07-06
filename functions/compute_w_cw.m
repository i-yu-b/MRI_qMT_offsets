function w_cw = compute_w_cw(pulse)
%{
w_cw = compute_w_cw(pulse) function computes mean omega1 value of shaped 
MT pulse using Sled and Pike CW approach: MT shaped pulse of duration t_mt 
approximated by cw rectangular pulse of TR duration and equivalent average 
power as shaped MT pulse
input:
    pulse - data structure containig all MT pulse properties
output:
    w_cw - averaged omega1 value
%}

TR = pulse.TR;
t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);
w_cw = sqrt(integral_value / TR);

end
