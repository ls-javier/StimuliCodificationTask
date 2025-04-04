function [p5, duration_stim] = pulses_5_mod(x,type_of_stim,amp)

l = length(x);

duration = floor(l/10);

if strcmp(type_of_stim, 'extended')
    time_pulse_2 = 2*duration;
    time_pulse_3 = 4*duration;
    time_pulse_4 = 6*duration;
    time_pulse_5 = 8*duration;

elseif strcmp(type_of_stim, 'grouped')
    time_pulse_2 = 3*duration;
    time_pulse_3 = 4*duration;
    time_pulse_4 = 5*duration;
    time_pulse_5 = 8*duration;
end

p5 = zeros(1,l); %definicion primer estimulo 

t = 0.5:duration-0.5;
p5(1,1:duration) =                           amp*sin(t*pi/duration); %shape of chain
p5(1,time_pulse_2+1:time_pulse_2+duration) = amp*sin(t*pi/duration);
p5(1,time_pulse_3+1:time_pulse_3+duration) = amp*sin(t*pi/duration);  
p5(1,time_pulse_4+1:time_pulse_4+duration) = amp*sin(t*pi/duration); 
p5(1,time_pulse_5+1:time_pulse_5+duration) = amp*sin(t*pi/duration);

duration_stim=[1,duration;time_pulse_2+1,time_pulse_2+duration;time_pulse_3+1,time_pulse_3+duration;...
    time_pulse_4+1,time_pulse_4+duration;time_pulse_5+1,time_pulse_5+duration]; % times of each pulses of the chain
