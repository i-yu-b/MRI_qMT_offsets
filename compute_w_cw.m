function w_cw = compute_w_cw(TR, pulse)
% Approximate MT pulse using Sled and Pike CW approach
% w_cw function computes cw rectangular pulse of TR duration and 
% equivalent average power as shaped MT pulse in pulse sequence 
% of duration t_mt

% input: TR - repetition time, pulse - set of pulse parameters
% output: averaged omega1 value

t_mt = pulse.t_mt;
omega1 = pulse.omega1;
integral_value = integral(omega1.^2, 0, t_mt);
w_cw = sqrt(integral_value / TR);

end


function w1rms = compute_w1rms(Pulse)
%compute_w1rms Compute the equivalent power of a rectangular pulse of same
%duration as the shaped pulse

Trf = Pulse.Trf;
omega2 = Pulse.omega2;
int = integral(omega2, 0, Trf);
w1rms = sqrt( int / Trf );

end