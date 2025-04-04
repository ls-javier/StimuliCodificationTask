function varargout = run_model_train(mode,T,task,vector_delay,histogram_ISI,varargin)
%Implementamos full-FORCE, para la red de espigas (LIF), la red auxiliar
%rate (v. continua) y el uso de RLS.

%cualquier modificación perdurará en todo el código
global N Ntilde n Nout Nin Nhint...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin u_tilde_fhint Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr ...
    randseed result_folder result_archivos_folder


%% Unpack varargin

uPC = varargin{1}; % varargin{n}: define n variables input
Imu = varargin{2};
hint = varargin{3};



%dimension of PC space
ntilde = size(uPC,2);        

if numel(varargin) > 3 %numel(): nº elementos del array
    dJ_PC = varargin{4};
    w = varargin{5};
    P = varargin{6};
else
% Si no hemos hecho el entrenamiento, tenemos que definir más variables.
    %learned spike net matrix
    dJ_PC = zeros(ntilde,2*n); %Equivalente al Ws de full-FORCE
    w = zeros(Nout,2*n); %output matrix % %Equivalente al W mayúscula de full-FORCE
    % n definido en linea 131 de demo_SN.m
    P = 1 * eye(2*n); %inverse covariance matrix; P(t) en RLS.             

end   
         
   


%% Figura 1 de la actividad de 5 neuronas junto con el esquema de los estímulos, targets...:

fh = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
'NumberTitle','off','ToolBar','none');

ah1 = axes(fh,'LineWidth',4,'FontSize',25,'Position',[0.1 0.71 0.8 0.2],...
'ylim',[-1.5 2.1]);
set(gca, 'FontWeight', 'bold')
ah2 = axes(fh,'LineWidth',4,'FontSize',25,'Position',[0.1 0.15 0.8 0.45]);
xlabel(ah2,'time [ms]');
ylabel(ah1,'Amplitude [mA]');
ylabel(ah2, 'Membrane Voltage [mV]');
set(gca, 'FontWeight', 'bold')

delete(findobj(ah1,'Type','Line'));
delete(findobj(ah2,'Type','Line'));

lh(1) = line(ah1,'Color','b','LineWidth',2,'Marker','none','LineStyle','-'); %output atractor
lh(2) = line(ah1,'Color','c','LineWidth',2,'Marker','none','LineStyle','-'); %output bump
lh(3) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 1
lh(4) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input central pulses
lh(5) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 5
lh(6) = line(ah1,'Color','r','LineWidth',2,'Marker','none','LineStyle','-'); %target atractor
lh(7) = line(ah1,'Color','m','LineWidth',2,'Marker','none','LineStyle','-'); %target bump

lh(8) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','--'); %pd
if strcmp(hint,'hint_decision')
    lh(9) = line(ah1,'Color','g','LineWidth',2,'Marker','none','LineStyle','--');
end              

nplt = 5;
for i = 1:nplt
    lhV(i) = line(ah2,'Color','k','LineWidth',2,'Marker','none','LineStyle','-');
end

        

%% Figura 2: RASTERS junto con el esquema de los estímulos, targets...:


fh2 = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
'NumberTitle','off','ToolBar','none');

ah3 = axes(fh2,'LineWidth',4,'FontSize',25,'Position',[0.1 0.71 0.8 0.2]);
set(gca, 'FontWeight', 'bold')
ah4 = axes(fh2,'LineWidth',4,'FontSize',25,'Position',[0.1 0.15 0.8 0.45]);
xlabel(ah4,'time [ms]');
ylabel(ah3,'Mean r(t) [mHz]');
ylabel(ah4,'Neurons');
set(gca, 'FontWeight', 'bold')

delete(findobj(ah3,'Type','Line'));
delete(findobj(ah4,'Type','Line'));

lh2(1) = line(ah3,'Color','b','LineWidth',2,'Marker','none','LineStyle','-'); 
lh2(2) = line(ah3,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 1
lh2(3) = line(ah3,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input central pulses
lh2(4) = line(ah3,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 5
lh2(5) = line(ah3,'Color','r','LineWidth',2,'Marker','none','LineStyle','-'); %target atractor
lh2(6) = line(ah3,'Color','m','LineWidth',2,'Marker','none','LineStyle','-'); %target bump

lh2(7) = line(ah3,'Color','k','LineWidth',2,'Marker','none','LineStyle','--'); %pd
if strcmp(hint,'hint_decision')
    lh2(8) = line(ah3,'Color','g','LineWidth',2,'Marker','none','LineStyle','--');
end       

lh_raster = line(ah4,'Color','k','LineWidth',1,'Marker','.','LineStyle','none');


%% Figura 3: Esquema de los estímulos, targets...:


fh3 = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
'NumberTitle','off','ToolBar','none');

ah5 = axes(fh3,'LineWidth',4,'FontSize',25,'Position',[0.1 0.15 0.8 0.8]);
xlabel(ah5,'time [ms]');
ylabel(ah5,'Amplitude [mA]');
set(gca, 'FontWeight', 'bold')

delete(findobj(ah5,'Type','Line'));

lh3(1) = line(ah5,'Color','b','LineWidth',2,'Marker','none','LineStyle','-'); %output atractor
lh3(2) = line(ah5,'Color','c','LineWidth',2,'Marker','none','LineStyle','-'); %output bump
lh3(3) = line(ah5,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 1
lh3(4) = line(ah5,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input central pulses
lh3(5) = line(ah5,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 5
lh3(6) = line(ah5,'Color','r','LineWidth',2,'Marker','none','LineStyle','-'); %target atractor
lh3(7) = line(ah5,'Color','m','LineWidth',2,'Marker','none','LineStyle','-'); %target bump

lh3(8) = line(ah5,'Color','k','LineWidth',2,'Marker','none','LineStyle','--'); %pd
if strcmp(hint,'hint_decision')
    lh3(9) = line(ah5,'Color','g','LineWidth',2,'Marker','none','LineStyle','--');
end 

   
%% Figura 4: RASTERS en periodo de equilibración:


fh4 = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
'NumberTitle','off','ToolBar','none');


ah6 = axes(fh4,'LineWidth',4,'FontSize',25,'Position',[0.1 0.15 0.8 0.8]);
xlabel(ah6,'time [ms]');
ylabel(ah6,'Neurons');
set(gca, 'FontWeight', 'bold')

delete(findobj(ah6,'Type','Line'));


lh_raster_Tinit = line(ah6,'Color','k','LineWidth',1,'Marker','.','LineStyle','none');


        

          
%% Initialize matrices for saving data

%nMSE, saved for each trial
nMSE = NaN(T,1);   % calcular el error (para cada ensayo) entre z_out y Ftilde 


EZ = 0; %used to compute running update of Error variance
NF = 0; %same, but target variance
  
       
% % % % % % %             %learned spike net matrix
% % % % % % %             dJ_PC = zeros(ntilde,2*n); %Equivalente al w minúscula de fullFORCE
% % % % % % %             w = zeros(Nout,2*n); %output matrix % %Equivalente al W mayúscula de fullFORCE
% % % % % % %             P = 1 * eye(2*n); %inverse covariance matrix    

%% Initalize state variables

%place network into the same state
%Fijamos la semilla para que la aleatoriedad de las variables sea siempre la misma
rng(randseed);  %% 3 es la semilla para el generador rng (en demo_SN por ej se define randseed=6, line 186)

x = 1e0 * randn(Ntilde,1); %teaching network state %target-generating network state equivalente a x_teach
zPC = zeros(ntilde,1); %state vector of learned dyanamics

%Variables que controlan el bucle (T -> Ensayos)
%time in current trial. Tiempo que transcurre para un ensayo. 
%Infinito para que ttrial>TTrial siempre (linea 218)
ttrial = inf;
%total time in current trial. Tiempo exacto que dura un ensayo (cte)
TTrial = 0;
%trial number. Contador de número de ensayos, define el final del bucle while go.
t_trial_num = 0;
%time across all trials
ttime = 1;
%flag to quit loop. Abre el bucle while (inicializador de ensayos)
go = true;


%product of PCs to Ntilde, Ntilde to N. For going directory from PCs to
%spiking network
uJ_uPC = uJ * uPC;  %% uJ: pesos u (matriz) en la ecuación (15) F_W(t). Transmite la info de rate aux a LIF principal.

%spike net state
V = 1e-3 * randn(N,1);

%learned input current
zJ = zeros(N,1);
%initial slow input current
z0s = zeros(N,1);
%inital fast input current
z0f = zeros(N,1);
%initial fast inhibition 
muglo = zeros(N,1);

% Son vectores de N componentes

%slow presynaptic current when using RLS
ss = zeros(n,1);
%fast presynaptic current when using RLS
sf = zeros(n,1);


taux=-dt/log(etaux); %No se usa
   

%% Run simulation

while go
    
    %generate trial data
    %Siempre ocurre porque ttrial=inf.
    if ttrial > TTrial   %% cuando el tiempo del trial (ttrial) supera la duracion TTrial del mismo
        
        %reset time for this trial to 1
        ttrial = 1;  
       
        %increase trial count by 1
        t_trial_num = t_trial_num + 1;
        
        %% Construimos el ensayo:
        % Define la duración y los estímulos
        [fin,fout,fhint,TTrial] = trial_conversion(vector_delay,hint,histogram_ISI);                   
        
        tspike = zeros(7*TTrial,1); %Storage variable for spike times. Indice de neuronas que ha disparado (qué neurona ha disparado)
        time_spike = zeros(7*TTrial,1); %Instante en el que una neurona dispara
        ns = 0; %Number of spikes, counts during simulation. Ver número de espigas.
         
         
        %Definimos ahora las conexiones que introducen algo en ambas redes.
        %fout input into each rate unit
        u_tilde_fout_fout = u_tilde_fout * fout;
        %fin input into each rate unit. 
        u_tilde_fin_fin = u_tilde_fin * fin;
        %fhint input into each rate unit. Pistas para la red, en función de cada tarea.
        u_tilde_fhint_fhint = u_tilde_fhint * fhint;

        %fin and fhint inputs into each spiking neuron. Pesos que entran en la red principal
        uJ_utilde_fin_fin = uJ * u_tilde_fin_fin;  
        
        if t_trial_num > Tinit   %% Tinit: el el nro de trial usadas para equilibrar
            
            %Make temporary empty arrays for saving things, or trial specific
            %inputs. Variables de almacenamiento
           
            %for collecting z(t), for plotting
            zs = zeros(Nout,TTrial);
                    
            %for collecting V(t), for plotting
            Vs = zeros(nplt,TTrial);

                    
                    
        end
        %Ya hemos construido el ensayo
    end
    %Ahora construimos el entrenamiento (importante)

    %rate model   (target-generating network)                                     
    fJ = Jtilde * tanh(x) + u_tilde_fout_fout(:,ttrial) + u_tilde_fhint_fhint(:,ttrial); %% ec6, P6
    xinf = fJ + u_tilde_fin_fin(:,ttrial);
    x = xinf + (x - xinf) * etaux;
    
    % x==tasa continua de disparo
    % 

    
    
    %spiking model             
    %project firing rate network inputs down to the PC basis
    fPC = uPC' * fJ; %%En 311: fJ --> fJs - En 403: temp = uPC' * fJs; 
    %Sólo metemos los input, no los targets
    Vinf = zJ + z0s + z0f - Imu ...   
        + muglo + uJ_utilde_fin_fin(:,ttrial);
    V = Vinf + (V - Vinf) * etaum;
        
    %ya hemos modificado el potencial de membrana
        
    index = find(V>=Vth);  %Find the neurons that have spiked 
    %definimos el indice de las neuronas que han disparado en este instante
    if ~isempty(index)     %Store spike times 
        tspike(ns+1:ns+length(index)) = index;
        time_spike(ns+1:ns+length(index))=0*index+ttrial;
        ns = ns + length(index);  % total number of spikes so far
    end  

    % Comando >=: 1 si ha disparado; 0 si no ha disparado
    S = V >= Vth;
    V(S) = Vth + Vr;                                  

    %Corrientes iniciales 
    z0s = z0s * etaus + sum(Js(:,S),2);
    z0f = z0f * etauf + sum(Jf(:,S),2);
    muglo = muglo * etauf + sum(Jmu(:,S),2); 

    %make presynaptic currents and concatenate

    % pregunta si las neuronas que han disparado tienen que ser modificadas (entrenadas)
    ss = etaus * ss + S(s_train);
    sf = etauf * sf + S(s_train);

    s = [ss;sf]; % s(t) en RLS (paper Parga, Serrano)
    %corriente que circula por las que vamos a modificar

        
    %Empezamos con el entrenamiento.
    %generate output and learned feedback output
    if t_trial_num > Tinit  % Muy importante


        %feedback currents in PC basis
        zPC = dJ_PC * s; % pesos que vamos a modificar por las corrientes modificables.
        %project from PC to spiking
        zJ = uJ_uPC * zPC; %proyeccion de las componentes principales. Cambio de base
        % para ponerlo en la misma dimensión que la red.

        %output
        z = w * s; %pesos output por la corriente que circula por esos pesos


    else

        %feedback is targets in PC basis. Define corriente como la que circula por la red principal
        zPC = fPC;
        %project from PC basis to spiking basis
        zJ = uJ_uPC * zPC;

        %output. Corriente z igual al target.
        z = fout(:,ttrial);

    end
        
    
    %after initial period
    if t_trial_num > Tinit
        
         %after initial period, save certain results for different computations
        
         %save for plotting and computing nMSE (mean squared error)
         %[]s -> variable de almacenamiento

         Vs(:,ttrial) = V(1:nplt);
         Vs(S(1:nplt),ttrial) = 2 * (Vth - Vr);

         zs(:,ttrial) = z;
                

         %Cálculo de la tasa instantanea de disparo de la población de neuronas:

        if ttrial == TTrial
            tasa=zeros(N,TTrial);
            for k=1:N
                ind=find(tspike==k);
                ind_time=time_spike(ind);
                tasa(k,ind_time)=1;
            end  

            % Cálculo de la tasa de disparo instantánea con la función firing rate.
            %primero calcula los tiempos de analisis, con una ventana de 250 ms que
            %se desplaza de a 50ms:
            t_ini=51;
            t_end=TTrial-t_ini;

            window_size = 50;

            slide=50;
            timeSamples = t_ini:slide:t_end; %window start time

            nTimes=length(timeSamples);

            [firing_rate_1] = firing_rate(tasa,timeSamples,window_size);
            mean=zeros(1,nTimes);
            for k=1:nTimes
                mean(1,k)=sum(firing_rate_1(:,k))/N;
            end  

            REC=zeros(1,TTrial);
                for k=1:nTimes
            REC(1,(k-1)*slide+1:(k-1)*slide+window_size+1)= mean(1,k);
            end

        end   


        %do RLS. Regla de aprendizaje
        if rand < 1/DTRLS % Definida en demo_SN, para controlar en qué instante entrenamos. (normalmente la mitad)

            xP = P * s;
            k = (1 + s' * xP)\xP';

            P = P - xP * k;
            dJ_PC = dJ_PC - (zPC - fPC) * k; % Modicamos pesos output con los pesos internos (eq. 34, Parga, Serrano)
            w = w - (z - fout(:,ttrial)) * k;

        end


        % text output, plot things, perform computations
        if ttrial == TTrial         %% AL FINAL DEL TRIAL


            %PERFORM COMPUTATIONS           

            %compute normalized output error on this trial
            nMSE(t_trial_num-Tinit) = sum(diag((zs - fout) * (zs - fout)'))/...
            sum(diag(fout * fout'));


            %printing and plotting
            %% Figura 1 de la actividad de 5 neuronas junto con el esquema de los estímulos, targets...:

            clc;
            fprintf('%s fullFORCE \nfullFORCE Error: %g\n%g trials of %g \n', ...
            mode, 100 * EZ/NF, t_trial_num-Tinit, T);

            set(lh(1),'XData',[ttrial-TTrial+1:ttrial],'YData',zs(1,:)); %output atractor
            set(lh(2),'XData',[ttrial-TTrial+1:ttrial],'YData',zs(2,:)); %output bump
            set(lh(3),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(1,:)); %input pulse 1
            set(lh(4),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(2,:)); %input central pulses
            set(lh(5),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(3,:)); %input pulse 5
            set(lh(6),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(1,:)); %target 1
            set(lh(7),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(2,:)); %target 2

            set(lh(8),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(4,:)); %pd
            if strcmp(hint,'hint_decision')
               set(lh(9),'XData',[ttrial-TTrial+1:ttrial],'YData',fhint(1,:));
            end    

            %plot some voltage trajectories         
            maxV = 0;
            for i = 1:nplt
                maxV = maxV + abs(min(Vs(i,:)));
                set(lhV(i),'XData',[ttrial-TTrial+1:ttrial],'YData', Vs(i,:) + maxV);
                maxV = maxV + max(Vs(i,:));
            end
            axis tight

            %% Figura 2: RASTERS junto con el esquema de los estímulos, targets...:
            clc;
            fprintf('%s fullFORCE \nfullFORCE Error: %g\n%g trials of %g \n', ...
            mode, 100 * EZ/NF, t_trial_num-Tinit, T);

            set(lh2(1),'XData',[ttrial-TTrial+1:ttrial],'YData',REC(1,:)*1000);
            set(lh2(2),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fin(1,:)); %input pulse 1
            set(lh2(3),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fin(2,:)); %input central pulses
            set(lh2(4),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fin(3,:)); %input pulse 5
            set(lh2(5),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fout(1,:)); %target 1
            set(lh2(6),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fout(2,:)); %target 2

            set(lh2(7),'XData',[ttrial-TTrial+1:ttrial],'YData',2*fin(4,:)); %pd
            if strcmp(hint,'hint_decision')
               set(lh2(8),'XData',[ttrial-TTrial+1:ttrial],'YData',fhint(1,:));
            end                             

            set(lh_raster,'XData',time_spike(:),'YData',tspike(:));


            %% Figura 3: Esquema de los estímulos, targets...:

            clc;
            fprintf('%s fullFORCE \nfullFORCE Error: %g\n%g trials of %g \n', ...
            mode, 100 * EZ/NF, t_trial_num-Tinit, T);

            set(lh3(1),'XData',[ttrial-TTrial+1:ttrial],'YData',zs(1,:)); %output atractor
            set(lh3(2),'XData',[ttrial-TTrial+1:ttrial],'YData',zs(2,:)); %output bump
            set(lh3(3),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(1,:)); %input pulse 1
            set(lh3(4),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(2,:)); %input central pulses
            set(lh3(5),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(3,:)); %input pulse 5
            set(lh3(6),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(1,:)); %target 1
            set(lh3(7),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(2,:)); %target 2

            set(lh3(8),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(4,:)); %pd
            if strcmp(hint,'hint_decision')
               set(lh3(9),'XData',[ttrial-TTrial+1:ttrial],'YData',fhint(1,:));
            end                          


            %% Drawnow and save figure                         

            drawnow;
            numero=t_trial_num-Tinit;

            cd([result_folder])

            if numero==2
                savefig(fh,'potncial_train_2')
                savefig(fh2,'raster_train_2')
                savefig(fh3,'esquema_train_2')
            elseif numero==50
                savefig(fh,'potncial_train_50')
                savefig(fh2,'raster_train_50')
                savefig(fh3,'esquema_train_50')
            elseif numero==100
                savefig(fh,'potncial_train_100')
                savefig(fh2,'raster_train_100')
                savefig(fh3,'esquema_train_100')
            elseif numero==150
                savefig(fh,'potncial_train_150')
                savefig(fh2,'raster_train_150')
                savefig(fh3,'esquema_train_150')
            elseif numero==200
                savefig(fh,'potncial_train_200')
                savefig(fh2,'raster_train_200')
                savefig(fh3,'esquema_train_200')
            elseif numero==250
                savefig(fh,'potncial_train_250')
                savefig(fh2,'raster_train_250')
                savefig(fh3,'esquema_train_250')
            end

            cd ..                        





        end

        %Definimos otros errores
        EZ = EZ + (z - fout(:,ttrial))' * (z - fout(:,ttrial));
        NF = NF + fout(:,ttrial)' * fout(:,ttrial);





        %counter of number of timesteps that have passed in total, over all
        %trials, only starts counting after initial period is passed
        ttime = ttime + 1;
        
    end
   
    
    if t_trial_num < Tinit
        if ttrial==TTrial

        set(lh_raster_Tinit,'XData',time_spike(:),'YData',tspike(:));

        drawnow;
        numero=t_trial_num;

        cd([result_folder])
        if numero==1
            savefig(fh4,'raster_train_Tinit_1')
        elseif numero==49
            savefig(fh4,'raster_train_Tinit_49')
        end   
        cd ..      
        end
        % Rasters = gráfica de puntos que define el nº de neuronas que han
        % disparado en cada instante.
    end
    
    %quit simulation loop
    if t_trial_num == T+Tinit && ttrial == TTrial 
        
        %quit loop
        go = false;                        %% <---------- fin del while go
        
    end
    
    %counter for number of timesteps that have passed in THIS trial
    ttrial = ttrial + 1;
    
end

%% output


varargout{1} = dJ_PC;
varargout{2} = w;
varargout{3} = P;

cd(result_archivos_folder)
    save('w','w')   
    save('dJ_PC','dJ_PC')  
    save('P','P')  
cd ..

delete(lh); %remove old line handles  
close(fh);
delete(lh2); %remove old line handles  
delete(lh_raster); %remove old line handles  
close(fh2);
delete(lh3); %remove old line handles 
close(fh3);
delete(lh_raster_Tinit); %remove old line handles  
close(fh4);

end
