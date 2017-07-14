function pulse_shape = sinc_pulse(t,t_mt, tbw)
%{
pulse_shape = sinc_pulse(t,t_mt,amp) function provides function handle to
describe lineshape(t) of sinc pulse of t_mt [s] duration and amplitude = 1

% inputs: 
%   t - time, s
%   t_mt - duration of pulse, s
% optional:
%   tbw - time-bandwidth product (total number of zero crossings, by
                                  default = 4)
% output: 
%   pulse_shape - function handle of computed pulse shape 
%}

if ~exist('tbw','var')
     % set default
     tbw = 4;
end
 
pulse_shape = sinc( tbw/t_mt * (t - t_mt/2) );
pulse_shape((t < 0 | t > t_mt)) = 0;

end
