function result = superlorentzian(delta,T2)
% superlorentzian lineshape function
% inputs: 
%   T2 - relaxation time of restiricted pool in s (scalar value)
%   delta - linear frequency offset of RF pulse, Hz (scalar value/vector)
% output: 
%   calculate lineshape (scalar value/vector)

integrand = @(u,delta,T2) sqrt(2/pi)*T2./...
              abs(3*u.^2-1).*exp(-2*((2*pi*delta*T2)'*(1./(3*u.^2-1))).^2);

result = integral(@(u)integrand(u,delta,T2),0,1,'ArrayValued',true)'; 

