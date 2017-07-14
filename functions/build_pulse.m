function pulse = build_pulse(metadata)
%{ 
function pulse = build_pulse(metadata) creates datastructure pulse, with
needed for fitting parameters of shaped MT pulse

input:
        dsfa
output:
        pulse - data structure with all MT shaped pulse parameters:
                pulse.t_mt - length of pulse, [s]
                pulse.tr - TR of pulse sequence, [s]
                pulse.amp - array of max pulse amplitudes, [Hz]    
                pulse.flip_angle - array of pulse flip angles, [degree]
                pulse.shape - string '<name>' describing shape of pulse
                currently available:
                  '<gauss>' - for gaussian pulse
                  '<sinc>' - for sinc pulse
                  '<hard>' - for hard pulse
                pulse.bw - bandwidth of pulse if available, [Hz]
                pulse.omega1_func - function handle describing shape of
                                    of pulse amplitude dependence on time 
                                    with max reached value = 1kHz
                pulse.omega1_sq_func - function handle describing shape of
                                    of pulse power dependence on time 
                                    with max reached value = (1kHz)^2
               
%}
    
% can be done only in case if one MT pulse per repetion is presented, i.e.
% metadata.PVM_MagTransPulsNumb ==1
    clear pulse;    
    gamma = 42.577478; % MHz/T
    
    % parcing all needed parameters of MT shaped pulse
    MTpulse_params = metadata.methpars.PVM_MagTransPulse1;
    pulse_params = strsplit(MTpulse_params,{', ','(',')'},...
                            'CollapseDelimiters',true); % cell string array 
    % length of MT pulse, [s]  
    pulse.t_mt = str2double(cell2mat(pulse_params(2)))/1000;

    % repetition time of pulse sequence
    pulse.TR = metadata.tr/1000; % [s]

    % amplitude of MT shaped pulse [Hz] 
    pulse.amp = [metadata.methpars.PVM_MagTransPower] * gamma;

    % flip angle of MT shaped pulse, [degree]
    pulse.flip_angle = str2double(cell2mat(pulse_params(4))); 

    % bandwidth of MT shaped pulse, [Hz]
    pulse.bw = str2double(cell2mat(pulse_params(3))); 

    % shape of MT pulse '<name of shape>'
    pulse.shape = metadata.methpars.PVM_MagTransPulse1Enum;

    %  setting shape function with standart amplitude = 1
    if strcmp(pulse.shape,'<gauss>') 
  
        pulse.omega1_func = @(t) gaussian_pulse(t, pulse.t_mt, pulse.bw);                   
        pulse.omega1_sq_func = @(t) gaussian_pulse(t, ...
                                                  pulse.t_mt, pulse.bw).^2;
    elseif strcmp(pulse.shape,'<sinc>') 
            % IMPORTANT! need to find time-bandwidth product for sinc 
            % pulse and set pulse.tbw = tbw and also pass tbw to
            % omega1_func. Right now it's set up to tbw = 4 by default
            pulse.omega1_func = @(t) sinc_pulse(t, pulse.t_mt);             
            pulse.omega1_sq_func = @(t) sinc_pulse(t, pulse.t_mt).^2;
        
    elseif strcmp(pulse.shape,'<hard>') 
            pulse.omega1_func = @(t) hard_pulse(t, pulse.t_mt);                          
            pulse.omega1_sq_func = @(t) hard_pulse(t, pulse.t_mt).^2;   
    end
end
    
    
    

    