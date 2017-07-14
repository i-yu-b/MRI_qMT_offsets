function pulse_shape = hard_pulse(t,t_mt)
%{
pulse_shape = hard_pulse(t,t_mt,amp) function provides function handle to
describe lineshape(t) of hard pulse of t_mt s duration

% inputs: 
%   t - time, s
%   t_mt - duration of pulse, s
%   amp - amplitude of pulse, s-1
% output: 
%   pulse_shape - function handle of computed pulse shape 
%}

pulse_shape = single(~(t < 0 | t > t_mt));
