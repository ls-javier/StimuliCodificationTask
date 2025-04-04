function varargout = run_model_test(mode,T,task,delay_value,histogram_ISI,varargin)

global N Ntilde n Nout Nin Nhint...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin u_tilde_fhint Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr ...
    randseed result_folder ...
    space_class


folder_result_delay = ['Test_Delay' num2str(delay_value)];

cd(result_folder)
    mkdir(folder_result_delay)
cd ..

%% Unpack varargin

uPC = varargin{1};
Imu = varargin{2};
dJ_PC = varargin{3};
w = varargin{4};
hint = varargin{5};


%% Figura 1 de la actividad de 5 neuronas junto con el esquema de los estï¿½mulos, targets...:

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
        
        

%% Figura 2: RASTERS junto con el esquema de los estï¿½mulos, targets...:


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

%% Figura 3: Esquema de los estï¿½mulos, targets...:


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


        
       
%% Figura 5: ISI histogram and ISI distribution:       

fh5 = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task %s: ISI histogram and ISI distribution. Period: %s ',task,histogram_ISI),...
'NumberTitle','off','ToolBar','none');

ah8 = axes(fh5,'LineWidth',4,'FontSize',18,'Position',[0.1 0.1 0.4 0.85]);
ah9 = axes(fh5,'LineWidth',4,'FontSize',18,'Position',[0.55 0.1 0.4 0.85]);
xlabel(ah8,'inter-spike interval [s]');
xlabel(ah9,'ISI [s]');
ylabel(ah8,'Count');
ylabel(ah9,'%');
set(gca, 'FontWeight', 'bold')
delete(findobj(ah8,'Type','Line'));
delete(findobj(ah9,'Type','Line'));

%         lh5 = line(ah8,'Color','b','LineWidth',1,'Marker','none','LineStyle','-');
lh6 = line(ah9,'Color','b','LineWidth',2,'Marker','none','LineStyle','-');




%% Initialize matrices for saving data


%nMSE, saved for each trial
nMSE = NaN(T,2);   

EZ = 0; %used to compute running update of Error variance
NF = 0; %same, but target variance

%% Initalize state variables

%place network into the same state
rng(randseed);  %% 3 es la semilla para el generador rng (en demo_SN por ej se define randseed=6, approx line 113)

ntilde = size(uPC,2);

x = 1e0 * randn(Ntilde,1); %teaching network state
zPC = zeros(ntilde,1); %state vector of learned dynamics

%time in current trial
ttrial = inf;
%total time in current trial
TTrial = 0;
%trial number
t_trial_num = 0;
%time across all trials
ttime = 1;
%flag to quit loop
go = true;

    
%product of PCs to Ntilde, Ntilde to N. For going directory from PCs to
%spiking network
uJ_uPC = uJ * uPC;  %% uJ

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


%slow presynaptic current when using RLS
ss = zeros(n,1);
%fast presynaptic current when using RLS
sf = zeros(n,1);

taux=-dt/log(etaux);


contador_errores=1;     %Contador usado para inicializar las variables de cálculo de trials correctos, línea 376 y 394.
%% Run simulation

while go
    
    %generate trial data
    if ttrial > TTrial   %% cuando el tiempo del trial (ttrial) supera la duracion TTrial del mismo
        
        %reset time for this trial to 1
        ttrial = 1;  
        
        %increase trial count by 1
        t_trial_num = t_trial_num + 1;
        

        %% Construimos el ensayos
        [fin,fout,fhint,TTrial,n_clase, delay_init, delay_fin,...
            t_tgt1_ini, t_tgt1_fin, t_tgt2_ini, t_tgt2_fin,...
            number_class, BAGAGE] = ...
            trial_conversion(delay_value, hint, histogram_ISI);

         if contador_errores==1
              contador_total=zeros(1,number_class);%Número total que se repite cada clase para el estandar en cada posicion
              test_correctos=zeros(1,number_class);%Número de tests correctos de cada clase.

             %Archivo usado para hacer el análisis PHASE SPACE.
             DATA_CONVERSION_1=cell(N/4,1);
             DATA_CONVERSION_2=cell(N/4,1);
             DATA_CONVERSION_3=cell(N/4,1);
             DATA_CONVERSION_4=cell(N/4,1);                 
             contador_trial=zeros(number_class,1);%número de clases

              contador_errores=0;
         end   
                                        
        
         
        tspike = zeros(7*TTrial,1); %Storage variable for spike times
        time_spike = zeros(7*TTrial,1);
        ns = 0; %Number of spikes, counts during simulation
         
        %fout input into each rate unit
        u_tilde_fout_fout = u_tilde_fout * fout;
        %fin input into each rate unit
        u_tilde_fin_fin = u_tilde_fin * fin;
        %fhint input into each rate unit
        u_tilde_fhint_fhint = u_tilde_fhint * fhint;
        
       
         %fin and fhint inputs into each spiking neuron
         uJ_utilde_fin_fin = uJ * u_tilde_fin_fin;  %% uJ
%          uJ_utilde_fhint_fhint = uJ * u_tilde_fhint_fhint;
         
        if t_trial_num > Tinit   %% Tinit: ¿el el nro de trial usadas para equilibrar?
            
            %Make temporary empty arrays for saving things, or trial specific
            %inputs
            
            %for collecting z(t), for plotting
            zs = zeros(Nout,TTrial);

            %for collecting V(t), for plotting
            Vs = zeros(nplt,TTrial);
                   
        end
    end
    
     %rate model   (target-generating network)                                          %% <----------
    fJ = Jtilde * tanh(x) + u_tilde_fout_fout(:,ttrial) + u_tilde_fhint_fhint(:,ttrial); %% ec6, P6
    xinf = fJ + u_tilde_fin_fin(:,ttrial);
    x = xinf + (x - xinf) * etaux;
    
    
    %spiking model                                          
    
    %project firing rate network inputs down to the PC basis
    fPC = uPC' * fJ; 

    Vinf = zJ + z0s + z0f - Imu ...  
           + muglo + uJ_utilde_fin_fin(:,ttrial);
    V = Vinf + (V - Vinf) * etaum;
        
        
    index = find(V>=Vth);  %Find the neurons that have spiked 
    if ~isempty(index)     %Store spike times (same as length(index)>0)
        tspike(ns+1:ns+length(index)) = index;
        time_spike(ns+1:ns+length(index))=0*index+ttrial;
        ns = ns + length(index);  % total number of spikes so far

        if ttrial==delay_init || ttrial==delay_init+1
            ns_delay_init=ns;

        elseif ttrial==delay_fin || ttrial==delay_fin+1

            ns_delay_fin=ns;

        end
    end  
        
    S = V >= Vth;
    V(S) = Vth + Vr;                                  

    z0s = z0s * etaus + sum(Js(:,S),2);
    z0f = z0f * etauf + sum(Jf(:,S),2);
    muglo = muglo * etauf + sum(Jmu(:,S),2);  

    %make presynaptic currents and concatenate


    ss = etaus * ss + S(s_train);
    sf = etauf * sf + S(s_train);

    s = [ss;sf];
          


    %generate output and learned feedback output
    if t_trial_num > Tinit

        %feedback currents in PC basis
        zPC = dJ_PC * s;
        %project from PC to spiking
        zJ = uJ_uPC * zPC;

        %output
        z = w * s;

    else

        %feedback is targets in PC basis
        zPC = fPC;
        %project from PC basis to spiking basis
        zJ = uJ_uPC * zPC;

        %output
        z = fout(:,ttrial);

    end
        
    %after initial period
    if t_trial_num > Tinit
        
        %after initial period, save certain results for different computations

        %save for plotting and computing nMSE
        Vs(:,ttrial) = V(1:nplt);
        Vs(S(1:nplt),ttrial) = 2 * (Vth - Vr);
        zs(:,ttrial) = z;
                
                
        % text output, plot things, perform computations
        if ttrial == TTrial         %% AL FINAL DEL TRIAL
            
            
            %% Cï¿½lculo de la tasa instantanea de disparo de la poblaciï¿½n de neuronas:

            tasa=zeros(N,TTrial);
            for k=1:N
               ind=find(tspike==k);
               ind_time=time_spike(ind);
               tasa(k,ind_time)=1;
            end  

            % Cï¿½lculo de la tasa de disparo instantï¿½nea con la funciï¿½n firing rate.
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

            
            
            
            %PERFORM COMPUTATIONS   
            %compute normalized output error on this trial
            nMSE(t_trial_num-Tinit) = sum(diag((zs - fout) * (zs - fout)'))/...
                                      sum(diag(fout * fout'));
            
            %% Para calcular el nï¿½mero y el porcentaje de test correctos:
            
            if Nout==1
                fout_values=fout;
                zs_values=zs;
            elseif Nout==2
                fout_values=fout(1,:);
                zs_values=zs(1,:);
            end
            %% Comprobar si el ensayo ha sido correcto para ambos targets (~434)
            %% Para el tgt 1
            fout1_cuadrado=fout_values(t_tgt1_ini:t_tgt1_fin).^2; 
            sum_fout1_cuadrado=sum(fout1_cuadrado);

            zs1_cuadrado=zs_values(t_tgt1_ini:t_tgt1_fin).^2;
            sum_zs1_cuadrado=sum(zs1_cuadrado);

            zs1_times_fout=zs_values.*fout_values;
            sumatorio_zs1_fout=sum(zs1_times_fout);

            correct_response1=sumatorio_zs1_fout/(sqrt(sum_zs1_cuadrado)*sqrt(sum_fout1_cuadrado)); 

            %% Para el tgt 2
            fout2_cuadrado=fout_values(t_tgt2_ini:t_tgt2_fin).^2; 
            sum_fout2_cuadrado=sum(fout2_cuadrado);

            zs2_cuadrado=zs_values(t_tgt2_ini:t_tgt2_fin).^2;
            sum_zs2_cuadrado=sum(zs2_cuadrado);

            zs2_times_fout=zs_values.*fout_values;
            sumatorio_zs2_fout=sum(zs2_times_fout);

            correct_response2=sumatorio_zs2_fout/(sqrt(sum_zs2_cuadrado)*sqrt(sum_fout2_cuadrado));
            
            HIT=0; %Si HIT=0, significa que la red da una respuesta incorrecta en este trial.
                   %Si HIT=1 significa que la red da la respuesta correcta.
                   
                         
            contador_total(n_clase)=contador_total(n_clase)+1;  

            if correct_response1 > 0.5 && correct_response2 > 0.5 %% d1 called different
                %%% Que el target sea una beta negativa (ngt=0) y que el ensayo sea
                %%% correcto significa que el output de la red se corresponde con el 
                %%% botón 'd1 called different'.

                HIT = 1;

                test_correctos(n_clase)=test_correctos(n_clase)+1;

            end                         


            %printing and plotting
            %% Figura 1 de la actividad de 5 neuronas junto con el esquema de los estï¿½mulos, targets...:
                        
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
                        
%% Figura 2: RASTERS junto con el esquema de los estï¿½mulos, targets...:
                        
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
                        

 
%% Figura 3: Esquema de los estï¿½mulos, targets...:

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
                        
      
%% Figura 5: ISI histogram and ISI distribution:   

            dim_index=length(tspike(ns_delay_init:ns_delay_fin));
            index_time=zeros(dim_index,2);

            index_time(:,1)=tspike(ns_delay_init:ns_delay_fin);
            index_time(:,2)=time_spike(ns_delay_init:ns_delay_fin);

%% ISI distribution     
            for jneuron=1:1:N                
                aux= find(index_time(:,1)==jneuron);
                both_spike_times = index_time(aux,2);
%                 spike_times = tspike(aux,:); % pone las dos columnas
                size1=length(both_spike_times);
                interspike_intervals = zeros(size1-1);
                interspike_intervals = both_spike_times(2:end) - both_spike_times(1:end-1);
%                 aux2=length(interspike_intervals)
                if jneuron==1
                   ISI_pop=interspike_intervals;
                else
                   ISI_pop=vertcat(ISI_pop,interspike_intervals); 
                end
            end
            
% plotting ISI distribution

            [dist, isi] = hist(ISI_pop);
            lh5=bar(ah8, isi, dist );


            dist=dist/sum(dist);
            set(lh6,'XData',isi,'YData',dist);
    

    
    
%% ARCHIVO PARA REALIZAR EL ANï¿½LISIS DE PHASE SPACE O INFORMACION MUTUA:

          
            for cl=1:number_class %tengo 4 clases.
                if n_clase==cl
                    counter = space_class * (cl - 1); %Si tengo 10 clases y hay un total de 1000 trials,
                                            % cada clase tendrá en torno a 100
                                            % trials. Pongo que cada clase se
                                            % espacie 150 índices de fila (con
                                            % el fin de ir almacenando).
                end
            end
            contador_trial(n_clase)=contador_trial(n_clase) + 1;
    
            for neuron=1:N

                 celda(neuron,contador_trial(n_clase) + counter,1)={n_clase};%Clase
                 celda(neuron,contador_trial(n_clase) + counter,2)={t_trial_num - Tinit};%Ensayo
                 celda(neuron,contador_trial(n_clase) + counter,3)={HIT};%HIT
                 celda(neuron,contador_trial(n_clase) + counter,4)={BAGAGE(7)};%P1 Estímulo posición 1
                 celda(neuron,contador_trial(n_clase) + counter,5)={BAGAGE(8)};%P2 Estímulo posición 2

                 celda(neuron,contador_trial(n_clase) + counter,7)={BAGAGE(1)};%PD
                 celda(neuron,contador_trial(n_clase) + counter,8)={BAGAGE(1)};%KD
                 celda(neuron,contador_trial(n_clase) + counter,9)={BAGAGE(2)};%SO1     
                 celda(neuron,contador_trial(n_clase) + counter,10)={BAGAGE(3)};%SO2
                 celda(neuron,contador_trial(n_clase) + counter,11)={BAGAGE(4)};%SF1
                 celda(neuron,contador_trial(n_clase) + counter,12)={BAGAGE(5)};%SF2    
                 celda(neuron,contador_trial(n_clase) + counter,13)={BAGAGE(6)};%KU
                 celda(neuron,contador_trial(n_clase) + counter,14)= {BAGAGE(6)};%PU
                 celda(neuron,contador_trial(n_clase) + counter,15)={0};%PK    
                 celda(neuron,contador_trial(n_clase) + counter,16)= {0};%RW
                 celda(neuron,contador_trial(n_clase) + counter,17)={0}; %A1
                 celda(neuron,contador_trial(n_clase) + counter,18)={0};%T1
                 celda(neuron,contador_trial(n_clase) + counter,19)= {0};%A2 
                 celda(neuron,contador_trial(n_clase) + counter,20)={0};%T2     
                 celda(neuron,contador_trial(n_clase) + counter,21)={BAGAGE(9)};%Inicio del target 
                 celda(neuron,contador_trial(n_clase) + counter,22)={BAGAGE(10)};%Fin del target 

                 celda(neuron,contador_trial(n_clase) + counter,23)={BAGAGE(11)};%Inicio del HINT de DECISIÓN 
                 celda(neuron,contador_trial(n_clase) + counter,24)={BAGAGE(12)};%Fin del HINT de DECISIÓN             

                 ind=find(tspike==neuron);
                 ind_time=time_spike(ind); 
                 celda(neuron,contador_trial(n_clase) + counter,6)={ind_time};%SPIKES TIMES 

            end             
        
            
%% Drawnow and save figure                         
                         
            drawnow;
            numero=t_trial_num-Tinit;
            cd([result_folder '/' folder_result_delay])

            if numero==2
                savefig(fh,['potncial_test_2_' num2str(delay_value)])
                savefig(fh2,['raster_test_2_' num2str(delay_value)])
                savefig(fh3,['esquema_test_2_' num2str(delay_value)])
                savefig(fh5,['histograma_test_2_' num2str(delay_value)])
            elseif numero==25
                savefig(fh,['potncial_test_25_' num2str(delay_value)])
                savefig(fh2,['raster_test_25_' num2str(delay_value)])
                savefig(fh3,['esquema_test_25_' num2str(delay_value)]) 
                savefig(fh5,['histograma_test_25_' num2str(delay_value)])
            elseif numero==50
                savefig(fh,['potncial_test_50_' num2str(delay_value)])
                savefig(fh2,['raster_test_50_' num2str(delay_value)])
                savefig(fh3,['esquema_test_50_' num2str(delay_value)])
                savefig(fh5,['histograma_test_50_' num2str(delay_value)])
            elseif numero==75
                savefig(fh,['potncial_test_75_' num2str(delay_value)])
                savefig(fh2,['raster_test_75_' num2str(delay_value)])
                savefig(fh3,['esquema_test_75_' num2str(delay_value)])
                savefig(fh5,['histograma_test_75_' num2str(delay_value)])
            elseif numero==100
                savefig(fh,['potncial_test_100_' num2str(delay_value)])
                savefig(fh2,['raster_test_100_' num2str(delay_value)])
                savefig(fh3,['esquema_test_100_' num2str(delay_value)])
                savefig(fh5,['histograma_test_100_' num2str(delay_value)])
            elseif numero==125
                savefig(fh,['potncial_test_125_' num2str(delay_value)])
                savefig(fh2,['raster_test_125_' num2str(delay_value)])
                savefig(fh3,['esquema_test_125_' num2str(delay_value)])
                savefig(fh5,['histograma_test_125_' num2str(delay_value)])
            elseif numero==150
                savefig(fh,['potncial_test_150_' num2str(delay_value)])
                savefig(fh2,['raster_test_150_' num2str(delay_value)])
                savefig(fh3,['esquema_test_150_' num2str(delay_value)])
                savefig(fh5,['histograma_test_150_' num2str(delay_value)])
            elseif numero==175
                savefig(fh,['potncial_test_175_' num2str(delay_value)])
                savefig(fh2,['raster_test_175_' num2str(delay_value)])
                savefig(fh3,['esquema_test_175_' num2str(delay_value)])
                savefig(fh5,['histograma_test_175_' num2str(delay_value)])
            elseif numero==200
                savefig(fh,['potncial_test_200_' num2str(delay_value)])
                savefig(fh2,['raster_test_200_' num2str(delay_value)])
                savefig(fh3,['esquema_test_200_' num2str(delay_value)])
                savefig(fh5,['histograma_test_200_' num2str(delay_value)])
            elseif numero==225
                savefig(fh,['potncial_test_225_' num2str(delay_value)])
                savefig(fh2,['raster_test_225_' num2str(delay_value)])
                savefig(fh3,['esquema_test_225_' num2str(delay_value)])
                savefig(fh5,['histograma_test_225_' num2str(delay_value)])
            elseif numero==250
                savefig(fh,['potncial_test_250_' num2str(delay_value)])
                savefig(fh2,['raster_test_250_' num2str(delay_value)])
                savefig(fh3,['esquema_test_250_' num2str(delay_value)])
                savefig(fh5,['histograma_test_250_' num2str(delay_value)])
            elseif numero==275
                savefig(fh,['potncial_test_275_' num2str(delay_value)])
                savefig(fh2,['raster_test_275_' num2str(delay_value)])
                savefig(fh3,['esquema_test_275_' num2str(delay_value)])
                savefig(fh5,['histograma_test_275_' num2str(delay_value)])
            end
            cd ../..
               
                          
        end
        

        EZ = EZ + (z - fout(:,ttrial))' * (z - fout(:,ttrial));
        NF = NF + fout(:,ttrial)' * fout(:,ttrial);

        %counter of number of timesteps that have passed in total, over all
        %trials, only starts counting after initial period is passed
        ttime = ttime + 1;
        
    end
    
    %quit simulation loop
    if t_trial_num == T+Tinit && ttrial == TTrial
        %quit loop
        go = false;                        %% <---------- fin del while go
    end
    
    %counter for number of timesteps that have passed in THIS trial
    ttrial = ttrial + 1;
    
end


%% Calcular porcentajes de clases

cd([result_folder '/' folder_result_delay])

diary('Accuracy.txt');

diary on;

cd ../..

porcentaje = test_correctos./contador_total.*100;

for i=1:number_class
disp(['Clase ' num2str(i)]) 
disp(['Número de ensayos correctos: ' num2str(test_correctos(i))])
disp(['Número total de ensayos: ' num2str(contador_total(i))])
disp(['Porcentaje de ensayos correctos: ' num2str(porcentaje(i))])
end 

diary off

%% Guardamos el archivo de datos para hacer el análisis de Phase Space y de MI:

suma_1=0;
suma_2=0;
for ge=1:number_class
suma_2=suma_2+contador_trial(ge);

celda_def(:, suma_1 + 1:suma_2, :) = celda(: , space_class * (ge - 1) + 1 : space_class * (ge - 1) + contador_trial(ge) , :);

suma_1=suma_1+contador_trial(ge);
end

celda=celda_def; %Redefinición.

[a1,b1,c1]=size(celda);

for k=1:a1
    Xo1=cell(b1+1,c1);
    Xo1(1,:)=[{'Clase'},{'N ensayo'},{'Hit'},{'P1'},{'P2'},{'spike times'},{'PD'},{'KD'},{'SO1'},{'SO2'},...
        {'SF1'},{'SF2'},{'KU'},{'PU'},{'PK'},{'RW'},{'A1'},{'T1'},{'A2'},{'T2'},...
        {'Target_ini'},{'Target_fin'},{'Hint_dec_ini'},{'Hint_dec_fin'}];
    for i=1:b1
        for j=1:c1 
            Xo1(i+1,j)=celda(k,i,j);
        end
    end

    if k<=N/4
        DATA_CONVERSION_1(k,1)={Xo1};
    elseif k>N/4 && k<=N/2
        DATA_CONVERSION_2(k-N/4,1)={Xo1};
    elseif k>N/2 && k<=3*N/4
        DATA_CONVERSION_3(k-2*N/4,1)={Xo1};  
    elseif k>3*N/4 
        DATA_CONVERSION_4(k-3*N/4,1)={Xo1};                 
    end
end

cd([result_folder '/' folder_result_delay])

save(['DATA_CONVERSION_' num2str(delay_value) '_1'],'DATA_CONVERSION_1')  
save(['DATA_CONVERSION_' num2str(delay_value) '_2'],'DATA_CONVERSION_2') 
save(['DATA_CONVERSION_' num2str(delay_value) '_3'],'DATA_CONVERSION_3')  
save(['DATA_CONVERSION_' num2str(delay_value) '_4'],'DATA_CONVERSION_4')   

clear DATA_CONVERSION_1; clear DATA_CONVERSION_2; clear DATA_CONVERSION_3; clear DATA_CONVERSION_4;

cd ../..
            

%% output

varargout{1} = nMSE;
delete(lh); %remove old line handles  
close(fh);
delete(lh2); %remove old line handles  
delete(lh_raster); %remove old line handles  
close(fh2);
delete(lh3); %remove old line handles 
close(fh3);
delete(lh5); %remove old line handles 
close(fh5);

end
