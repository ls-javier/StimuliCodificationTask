function varargout = run_model_PCA(mode,T,task,hint,varargin)

global N Ntilde n Nout Nin Nhint...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin u_tilde_fhint Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr

%% bucle temporal (while go): linea 166; calcula fJ: 228; 
%% calculo de las PCs (uPC): 356 ; projecta input fJ sobre PC: 236, tb: 403
%% din de x (rate): 231; din de V (spikes): 241

%% RLS: 320

%% Figure

fh = figure('Color','w', 'menubar', 'none', 'Name', sprintf('task: %s',task),...
    'NumberTitle','off','ToolBar','none');

ah1 = axes(fh,'LineWidth',2,'FontSize',12,'Position',[0.1 0.7 0.8 0.2],...
    'ylim',[-1.1 2.1]);
ah2 = axes(fh,'LineWidth',2,'FontSize',12,'Position',[0.1 0.1 0.8 0.5]);
xlabel(ah2,'time (s)');

delete(findobj(ah1,'Type','Line'));
delete(findobj(ah2,'Type','Line'));

lh(1) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 1
lh(2) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input central pulses
lh(3) = line(ah1,'Color','k','LineWidth',2,'Marker','none','LineStyle','-'); %input pulse 5
lh(4) = line(ah1,'Color','r','LineWidth',2,'Marker','none','LineStyle','-'); %target atractor
lh(5) = line(ah1,'Color','m','LineWidth',2,'Marker','none','LineStyle','-'); %target bump
lh(6) = line(ah1,'Color','g','LineWidth',1,'Marker','none','LineStyle','--'); %hint

    nplt = 5;
    for i = 1:nplt
        lhr(i) = line(ah2,'Color','r','LineWidth',2,'Marker','none','LineStyle','-');
    end

%% Initialize matrices for saving data

    %for saving covariance of firing rate network inputs
    COV_fJ = zeros(Ntilde);
    


%% Initalize state variables

%place network into the same state
rng(3);  %% 3 es la semilla para el generador rng (en demo_SN por ej se define randseed=6, approx line 113)

x = 1e0 * randn(Ntilde,1); %teaching network state %target-generating network state equivalente a x_teach
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


taux=-dt/log(etaux);

%% Run simulation

while go
    
    %generate trial data
    if ttrial > TTrial   %% cuando el tiempo del trial (ttrial) supera la duracion TTrial del mismo
        
        %reset time for this trial to 1
        ttrial = 1;  
        
        %increase trial count by 1
        t_trial_num = t_trial_num + 1;
        
       switch task
            case 'osc4'
                [fin,fout,TTrial,fhint] = trial_osc4(dt);   
            case 'ready_set_go'
                [fin,fout,fhint,TTrial] = trial_ready_set_go(dt,taux,hint);
       end
        
        
        %fout input into each rate unit
        u_tilde_fout_fout = u_tilde_fout * fout;
        %fin input into each rate unit
        u_tilde_fin_fin = u_tilde_fin * fin;
         %fhint input into each rate unit
        u_tilde_fhint_fhint = u_tilde_fhint * fhint;
        
      
        if t_trial_num > Tinit   %% Tinit: ¿el el nro de trial usadas para equilibrar?
            
            %Make temporary empty arrays for saving things, or trial specific
            %inputs
           
                    %for collecting fJ of the rate network, for doing PCA  
                    fJs = zeros(Ntilde,TTrial);             %% <----------
                    
               
            
        end
        
    end
    
    
      %rate model   (target-generating network)        
      fJ = Jtilde * tanh(x) + u_tilde_fout_fout(:,ttrial) + u_tilde_fhint_fhint(:,ttrial); %% ec6, P6
    xinf = fJ + u_tilde_fin_fin(:,ttrial);
    x = xinf + (x - xinf) * etaux;
    
 
    
    
    %after initial period
    if t_trial_num > Tinit
        
        %after initial period, save certain results for different computations
       
                %save for plotting and computing cov matrix
                fJs(:,ttrial) = fJ;
                
            
        
        % text output, plot things, perform computations
        if ttrial == TTrial         %% AL FINAL DEL TRIAL?
            
            %PERFORM COMPUTATIONS            %% <---------- %% <----------
                                  
                    
                    %update covariance matrix
                    COV_fJ = COV_fJ + fJs * fJs';  %% En 311: fJ --> fJs
                    
                    %compute PCs from cov matrix, and arrange them from
                    %largest to smallest
                    [uPC,d] = eig(COV_fJ);
                    uPC = fliplr(uPC);
                    d = flipud(diag(d));
                    
                    %pick the number of dimensions that accounts 
                    %for fracPCA of variance            %% <----------
                    ntilde = find(cumsum(d)/sum(d) > fracPCA,1, 'first');
                    uPC = uPC(:,1:ntilde);
                    
               
            
            %printing and plotting
            clc;

            %plot fout and fin
            set(lh(1),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(1,:)); %input pulse 1
            set(lh(2),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(2,:)); %input central pulses
            set(lh(3),'XData',[ttrial-TTrial+1:ttrial],'YData',fin(3,:)); %input pulse 5
            set(lh(4),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(1,:)); %target 1
            set(lh(5),'XData',[ttrial-TTrial+1:ttrial],'YData',fout(2,:)); %target 2
            set(lh(6),'XData',dt*[ttrial-TTrial+1:ttrial],'YData',fhint(1,:));
                    
                    fprintf('%s, %g PCs, %g trials of %g \n', mode, ntilde, t_trial_num-Tinit, T);
                    
                    %project down to the PCs
                    temp = uPC' * fJs;
                    max_temp = 0;
                    %plot a few of them
                    for i = 1:nplt
                        max_temp = max_temp + abs(min(temp(i,:)));
                        set(lhr(i),'XData',dt*[ttrial-TTrial+1:ttrial],'YData', temp(i,:) + max_temp);
                        max_temp = max_temp + max(temp(i,:));
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
        go = false;                        %% <---------- fin del while go
        
    end
    
    %counter for number of timesteps that have passed in THIS trial
    ttrial = ttrial + 1;
    
end

%% output
        varargout{1} = uPC;
 
close(fh);
