function W = calculate_W(delta,T2r,lineshape)
% calculate saturation rate
switch lineshape
    case 'gaussian'
        W = gaussian(delta, T2r);
    case 'lorentzian'
        W = lorentzian(delta, T2r);
    case 'superlorentzian'
        W = superlorentzian(delta, T2r);
    case 'superlorentzian_res'
        W = superlorentzian_res(delta, T2r);
end

end


function W = computeW(G, Pulse)
% Compute Mean saturation rate <W(delta,alpha)> for given G(delta)
% ----------------------------------------------------------------------------------------------------
% Written by: Jean-Fran?ois Cabana, 2016
% ----------------------------------------------------------------------------------------------------


W = zeros(length(Pulse),1);

for ii = 1:length(Pulse)
    omega2 = Pulse(ii).omega2;
    Trf = Pulse(ii).Trf;
    int = integral(omega2, 0, Trf, 'ArrayValued', true);
    W(ii) = G * pi/Trf * int;
end