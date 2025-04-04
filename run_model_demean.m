function varargout = run_model_demean(mode,T,task,vector_delay,histogram_ISI,varargin)

global N Ntilde n Nout Nin Nhint...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin u_tilde_fhint Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr ...
    randseed result_folder result_archivos_folder


%% Unpack varargin:

uPC = varargin{1};
hint = varargin{2};
    
  


%% Figura 1 de la actividad de 5 neuronas junto con el esquema de los estímulos, targets...:

    fh = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
    'NumberTitle','off','ToolBar','none');

    ah1 = axes(fh,'LineWidth',4,'FontSize',25,'Position',[0.1 0.71 0.8 0.2],...
    'ylim',[-1.5 2.1]);
    set(gca, 'FontWeight', 'bold')
    ah2 = axes(fh,'LineWidth',4,'FontSize',25,'Position',[0.1 0.15 0.8 0.45]);
    xlabel(ah2,'time [ms]');
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
    
    
       


%% Initialize matrices for saving data

    
    %computed mean of input into each spiking neuron due to learned signal
    sum_Imu = zeros(N,1);
    Imu = zeros(N,1);
    
    %mean firing rate, saved for each trial (in demean mode)
    mFR = NaN(T,1);
    
    EZ = 0; %used to compute running update of Error variance
    NF = 0; %same, but target variance
        
    %dimension of PC space
    ntilde = size(uPC,2);
%% Initalize state variables

%place network into the same state
rng(randseed);  %% 3 es la semilla para el generador rng (en demo_SN por ej se define randseed=6, approx line 113)

x = 1e0 * randn(Ntilde,1); %teaching network state %target-generating network state equivalente a x_teach
zPC = zeros(ntilde,1); %state vector of learned dyanamics

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
    uJ_uPC = uJ * uPC;  %% uJ: ¿es la u en ec 6, P6?
    
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
        

    taux=-dt/log(etaux);
%% Run simulation

while go
    
    %generate trial data
    if ttrial > TTrial   %% cuando el tiempo del trial (ttrial) supera la duracion TTrial del mismo
        
        %reset time for this trial to 1
        ttrial = 1;  
        
        %increase trial count by 1
        t_trial_num = t_trial_num + 1;
        
        %% Construimos el ensayo:
        [fin,fout,fhint,TTrial] = trial_conversion(vector_delay,hint,histogram_ISI);             
        
        %fout input into each rate unit
        u_tilde_fout_fout = u_tilde_fout * fout;
        %fin input into each rate unit
        u_tilde_fin_fin = u_tilde_fin * fin;
        %fhint input into each rate unit
        u_tilde_fhint_fhint = u_tilde_fhint * fhint;
        
       
         %fin and fhint inputs into each spiking neuron
         uJ_utilde_fin_fin = uJ * u_tilde_fin_fin;  %% uJ: ¿es la u en ec 6, P6?
%          uJ_utilde_fhint_fhint = uJ * u_tilde_fhint_fhint;
        
        if t_trial_num > Tinit   %% Tinit: ¿el el nro de trial usadas para equilibrar?
            
            %Make temporary empty arrays for saving things, or trial specific
            %inputs
            
                    
            %for collecting z(t), for plotting
            zs = zeros(Nout,TTrial);

            %for collecting V(t), for plotting
            Vs = zeros(nplt,TTrial);
              
                    
            %for collecting zJ,z0f and z0s, for mean
            %computation
            all_zs = zeros(N,TTrial);
            %for collecting the spikes, for computing the mean
            %firing rate
            Ss = zeros(N,TTrial);
                       
                 
            
        end
        
    end
    
    %rate model   (target-generating network)                                          %% <----------
    fJ = Jtilde * tanh(x) + u_tilde_fout_fout(:,ttrial) + u_tilde_fhint_fhint(:,ttrial); %% ec6, P6
    xinf = fJ + u_tilde_fin_fin(:,ttrial);
    x = xinf + (x - xinf) * etaux;
    
    
    %spiking model                      
    %project firing rate network inputs down to the PC basis
    fPC = uPC' * fJ; %%En 311: fJ --> fJs - En 403: temp = uPC' * fJs; 


    Vinf = zJ + z0s + z0f - Imu ...   
           + muglo + uJ_utilde_fin_fin(:,ttrial);
    V = Vinf + (V - Vinf) * etaum;

    S = V >= Vth;
    V(S) = Vth + Vr;                             

    z0s = z0s * etaus + sum(Js(:,S),2);
    z0f = z0f * etauf + sum(Jf(:,S),2);
    muglo = muglo * etauf + sum(Jmu(:,S),2);  
        
        
    %generate output and learned feedback output
    if t_trial_num > Tinit

        %use target as feedback 
        zPC = fPC;
        %project from PC to spiking basis
        zJ = uJ_uPC * fPC;

        %output is just the target
        z = fout(:,ttrial);



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
                

        %save for mean computation, project targets from PC
        %basis to spiking basis, and add initial slow and fast
        %currents
        all_zs(:,ttrial) = uJ_uPC * fPC + z0s + z0f;
        %save for mean firing rate computation
        Ss(:,ttrial) = S;
                       
            
           
        
        % text output, plot things, perform computations
        if ttrial == TTrial      
            
            %PERFORM COMPUTATIONS       
           

            %add sum of current trial's currents to running sum
            sum_Imu = sum_Imu + sum(all_zs,2); %% all_zs: 339
            %divde by total elapsed time to compute mean
            Imu = sum_Imu/ttime;

            %sum all spikes on this trial, and normalize by number
            %of neurons and elapsed time on this trial to compute
            %mean population firing rate on this trial
            mFR(t_trial_num-Tinit) = sum(Ss(:))/(N*TTrial*dt);
                  
            
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
             
            
            %plot generated output

            fprintf('%s, %g trials of %g \n', mode, t_trial_num-Tinit, T);
                        

            %plot some voltage trajectories      
            maxV = 0;
            for i = 1:nplt
                maxV = maxV + abs(min(Vs(i,:)));
                set(lhV(i),'XData',[ttrial-TTrial+1:ttrial],'YData', Vs(i,:) + maxV);
                maxV = maxV + max(Vs(i,:));
            end
            axis tight
                    
            drawnow;
                          
        end
        
        %counter of number of timesteps that have passed in total, over all
        %trials, only starts counting after initial period is passed
        ttime = ttime + 1;
        
    end
    
    %quit simulation loop
    if t_trial_num == T+Tinit && ttrial == TTrial
        
        %quit loop
        go = false;  %% <---------- fin del while go
        
    end
    
    %counter for number of timesteps that have passed in THIS trial
    ttrial = ttrial + 1;
    
end

%% output

varargout{1} = Imu;
varargout{2} = mFR;

cd(result_archivos_folder)
	save('Imu','Imu')
cd ..   
        
close(fh);
