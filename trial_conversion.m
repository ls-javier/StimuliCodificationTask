function [fin, fout, fhint, TTrial, n_clase, delay_init, delay_fin, ...
          t_tgt1_ini, t_tgt1_fin, t_tgt2_ini, t_tgt2_fin, ...
          number_class, BAGAGE] = ...
    trial_conversion(vector_delay, hint, histogram_ISI)

global Nout Nin Nhint

%% Inicializamos las variables necesarias

number_class = 4;

pd = 40; % Duración del Probe Down.
pu = 40; % Duración del Probe Up.

delay = 1500; % Time bt fin_stim_2 and pu
pre_stim1=1000;% Time to first stim.
stim_l = 1000; % Time for stim
inter_stim=1500; %Time bt stimuli
decision_l=1000; % Time for decission.
inter_trial=1000;% Time between trials.

tgt_amp = 1; % Amplitude of target.
stim_amp = 1; % Amplitude of stimulus

%% Construimos las clases, pu y pd

if length(vector_delay) == 1
    inter_stim = vector_delay;
else
    inter_stim = vector_delay(randi(length(vector_delay)));
    disp(['Delay ' num2str(inter_stim)])
end

n_clase=randi([1,4]); % Elige al azar la clase.

TTrial = pd + pre_stim1 + inter_stim + 2*stim_l + inter_trial;
t_stim1_ini= pd + pre_stim1;
t_stim1_fin=t_stim1_ini+stim_l;
        
t_stim2_ini=t_stim1_fin+inter_stim;
t_stim2_fin=t_stim2_ini+stim_l;

ninp1=[1,1,0,0];% MMM indica si en la clase n el patron 1 se aplica (1) o no (0) durante el primer periodo de estimulacion 
ninp2=[1,0,1,0];% MMM indica si en la clase n el patron 1 se aplica (1) o no (0) durante el 2do periodo de estimulacion  
                % MMM El cero se tomara como presencia del patron 2 (V mas abajo la construccion de los inputs 1 y 2)

% ninp1(n)indica si en la clase n el patron "long" (Grouped, L) se aplica (1) o no (0) durante la 1er estimulacion
% ninp2(n)indica si en la clase n el patron "long" (Grouped, L) se aplica (1) o no (0) durante la 2da estimulacion


%% Construcción del input de la función

fin=zeros(Nin,TTrial); %fila 1= estimulos input; Fila 2= pu pd. Dim= 2 input x duracion ensayo

%% Primer estímulo
type_of_stim1 = 'extended';
if ninp1(n_clase)==1
    type_of_stim1 = 'grouped';
end

[p5,duration_stim]=pulses_5_mod(fin(1,t_stim1_ini:t_stim1_fin),type_of_stim1,stim_amp); %Definimos cadena de pulsos para 1er estimulo

duration_stim_mod=duration_stim+t_stim1_ini;

fin(1,duration_stim_mod(1,1):duration_stim_mod(1,2))=p5(1,duration_stim(1,1):duration_stim(1,2));
fin(2,duration_stim_mod(2,1):duration_stim_mod(4,2))=p5(1,duration_stim(2,1):duration_stim(4,2));
fin(3,duration_stim_mod(5,1):duration_stim_mod(5,2))=p5(1,duration_stim(5,1):duration_stim(5,2));


%% Redefinición de fin para el segundo estímulo
type_of_stim2 = 'extended';
if ninp2(n_clase)==1
    type_of_stim2 = 'grouped'; 
end

[p5,duration_stim]=pulses_5_mod(fin(1,t_stim2_ini:t_stim2_fin),type_of_stim2,stim_amp); %Definimos cadena de pulsos para 1er estimulo

duration_stim_mod=duration_stim+t_stim2_ini;

fin(1,duration_stim_mod(1,1):duration_stim_mod(1,2))=p5(1,duration_stim(1,1):duration_stim(1,2));
fin(2,duration_stim_mod(2,1):duration_stim_mod(4,2))=p5(1,duration_stim(2,1):duration_stim(4,2));
fin(3,duration_stim_mod(5,1):duration_stim_mod(5,2))=p5(1,duration_stim(5,1):duration_stim(5,2));

fin(4,1:pu)=stim_amp; %Es el PD.

%% 1. Construimos el primer target (atractor)      

fout=zeros(Nout,TTrial);

t_tgt1_ini=t_stim1_ini;
[~,t_tgt1_fin, ~, ~]=pulses_5(fin(1,t_stim1_ini:t_stim1_fin),type_of_stim1,stim_amp);

t_tgt2_ini=t_stim2_ini;
[~,t_tgt2_fin,~, ~]=pulses_5(fin(1,t_stim2_ini:t_stim2_fin),type_of_stim2,stim_amp);


fout(1,t_tgt1_ini:t_tgt1_ini+t_tgt1_fin)=1; % MMM el tgt durante la primer estimulacion
fout(1,t_tgt2_ini:t_tgt2_ini+t_tgt2_fin)=1;

%% 2. Construimos el segundo target (bump)

[~, ~, t_tgt3_ini, t_tgt3_fin] = pulses_5(fin(1,t_stim1_ini:t_stim1_fin),type_of_stim1,stim_amp); % inicio del estimulo

t_tgt3_inicial=t_tgt3_ini+t_stim1_ini;  %Cambiamos de sistema de referencia al inicio del ensayo
t_tgt3_final=t_tgt3_inicial+t_tgt3_fin;

[~, ~, t_tgt4_ini, t_tgt4_fin] = pulses_5(fin(1,t_stim2_ini:t_stim2_fin),type_of_stim2,stim_amp);  % final del estimulo

t_tgt4_inicial=t_stim2_ini+t_tgt4_ini;
t_tgt4_final=t_tgt4_inicial + t_tgt4_fin;


longitud_1=length(t_tgt3_inicial:t_tgt3_final);
longitud_2=length(t_tgt4_inicial:t_tgt4_final);  

if ninp1(n_clase) == 1  % solo habrá bump en el primer estimulo si es agrupado
    IO1 = tgt_amp* betapdf(linspace(0,1,longitud_1),4,4)/max(betapdf(linspace(0,1,longitud_1),4,4));
else
    IO1=0;
end

if ninp2(n_clase) == 1  % lo mismo para el segundo estímulo
    IO2 = tgt_amp* betapdf(linspace(0,1,longitud_2),4,4)/max(betapdf(linspace(0,1,longitud_2),4,4));
else
    IO2=0;
end

fout(2,t_tgt3_inicial:t_tgt3_final)=IO1;
fout(2,t_tgt4_inicial:t_tgt4_final)=IO2;

%% Construimos la función fhint:
if strcmp(hint,'nohint')   % De momento estamos en este caso
   fhint = zeros(Nhint,TTrial);
   t_hint_decision_ini = NaN;
   t_hint_decision_fin = NaN;
elseif strcmp(hint,'hint_decision') 
   fhint = zeros(Nhint,TTrial);
   % HINT DE DECISIÓN DESPUES DEL SEGUNDO ESTÍMULO:
   
   t_hint_decision_ini=t_stim2_fin + 101;
   t_hint_decision_fin=t_stim2_fin + 101 + decision_l;

% %    t_hint_decision_ini=t_stim2_ini + 445 + 1;
% %    t_hint_decision_fin=t_stim2_ini + 445 + 1 + decision_l;

% %    t_hint_decision_ini= t_stim2_ini + estandar + 1;
% %    t_hint_decision_fin= t_hint_decision_ini + decision_l;
   
   IO3_hint = ((2*ngt-1)*tgt_amp)* betapdf(linspace(0,1,decision_l+1),4,4)/max(betapdf(linspace(0,1,decision_l+1),4,4));%target output
   
   fhint(1,t_hint_decision_ini:t_hint_decision_fin)=IO3_hint; 

end



%% Variables necesarias para construir el Histograma:
switch histogram_ISI
       case 'target'
            delay_init=t_tgt1_ini-1;
            delay_fin=t_tgt1_fin-1;
            
       case 'delay'
            delay_init=t_stim1_fin-1;
            delay_fin=t_stim2_ini-1;            
end

 fprintf('Estimulo en posicion 1, P1: %g ms \nEstimulo en posicion 2, P2: %g ms\n', ninp1,ninp2);   

 BAGAGE = [0,t_stim1_ini,t_stim2_ini,t_stim1_fin,t_stim2_fin,t_stim2_fin+delay,ninp1,ninp2,...
           t_tgt1_ini, t_tgt1_fin, t_tgt2_ini, t_tgt2_fin,t_hint_decision_ini,t_hint_decision_fin];

end

