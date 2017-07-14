function dMf = Bloch_eq(t, Mf, T2f, pulse_shape, offset, amp)
%{
 function function dMf = Bloch_eq(t, Mf, T2f, pulse, offsets)
 calculates Bloch differential equations in case of  RF excitation 
 and no MT effect
input:
        Mf - free pool magnetization vector Mf = [Mxf, Myf, Mzf]
        pulse_shape - function handle describing shape of RF pulse (t)
        amp - amplitude of pulse
        T2f - relaxation time of free pool, [s-1]
        offset - frequency offsets for RF pulse, [Hz]
output:
        dMf - delta value of free pool magnetization vector [Mxf, Myf, Mzf]
%}
omega1 = pulse_shape(t) * amp;
dMf    =  zeros(size(Mf));
dMf(1) = - Mf(1)/T2f - 2 * pi * offset * Mf(2);
dMf(2) = - Mf(2)/T2f + 2 *pi * offset * Mf(1) + omega1 * Mf(3);
dMf(3) = - omega1 * Mf(2);
