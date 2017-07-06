function w_rp = compute_w_rp(pulse)
%{
w_rp = compute_w_rp(pulse) function computes mean omega1 value of shaped 
MT pulse using Sled and Pike RP approach: MT shaped pulse of duration t_mt 
approximated by cw rectangular pulse of duration equals to 
full-width-at-half-maximum of the shaped MT pulse and equivalent average 
power

input:
    pulse - data structure containig all MT pulse properties
output:
    w_rp - averaged omega1 value
%}

TR = pulse.TR;
t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);

if strcmp(pulse.shape,'hard')
    tau = t_mt;
else
    t = 0:2000:t_mt;
    tau = fwhm(t,omega1);
end

w_rp = sqrt(integral_value/tau);

end
