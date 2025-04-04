clc
clear
%This script will train a network of LIF neurons based on the method
%described in "Using firing rate dynamics to train recurrent networks of
%spiking model neurons" (https://arxiv.org/abs/1601.07620)

%1. PCs of the firing rate dyanmics are computed, so that only the most
%dominate signals are learned in the spiking network.
%2. The mean input into each spiking neuron is computed and subtracted, to
%control aphysiological firing rates.
%3. RLS is used to train the inputs into each spiking neuron.
%4. Training performance is tested.

global N Ntilde n Nout Nin Nhint...
    dt Tinit ...
    fracPCA DTRLS ...
    s_train u_tilde_fout u_tilde_fin u_tilde_fhint Jtilde uJ Jmu Jf Js ...
    etaux etauf etaus etaum Vth Vr ...
    randseed ...
    result_folder result_archivos_folder ...
    space_class



%% ¿HACER LA ETAPA TRAIN? ¿O PASAR DIRECTAMENTE AL TEST DESCARGANDO LOS ARCHIVOS NECESARIOS?
doing_train = 'YES';
% % doing_train = 'NO';


if any(strcmp(doing_train,{'YES'}))

    result_folder = 'Resultados'; %Directorio donde se encuentran los archivos de datos de todas las redes.
    result_archivos_folder = 'Archivos';
    
    mkdir(result_folder)

    mkdir(result_archivos_folder)
    
elseif any(strcmp(doing_train,{'NO'}))
    result_folder = 'Result_wout_train'; %Directorio donde se encuentran los archivos de datos de todas las redes.
    result_archivos_folder = 'Archivos'; % Este archivo ya esta creado y contiene los diferentes arrays
    
    mkdir(result_folder)

end


%% Number of neurons

%N: number of spiking neurons
%Ntilde: number of rate units

N = 500;
Ntilde = 500;

%% How long to run each routine and which ones to run

%init: number of trials for dynamics to settle from initial condition
%PCA: how many trials to compute PCs
%demean: how many trials to compute mean
%RLS: how many trials to do RLS training for
%test: how many trials to test for, after training

Tinit = 50; % Necesario para que el sistema evolucione hacia un estado más estable
% % % TPCA = 300;
Tdemean = 1000; %1000 

TRLS = 4000; %1000 para comprobar (o un poco menos), 5000 para un buen entrenamiento

Ttest = 4000; %100 para comprobar, 4000 para un buen test

%% Valor del delay entre d1 y d2:
delay_values = 1500;
% % % delay_values = [1000,1500,2000,2500];

%% Gain parameters

%gtilde: gain of recurrent connectivity
%gfout: gain of output feedback
%gfin: gain of input signal
%gz: learned feedback gain of spiking net
%muf: global inhibition
%gs: gain of slow random synapses
%gf: gain of fast random synpases

gtilde = 1.4; %También es el gr en fullFORCE
gfout = 1.0;
gfin = 1.0;
gfhint=1.0;
gz = 6;
muf = -0.3;

% gs = 0.11; 
% gf = 0.13;

gs = 0;
gf = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% task
%%%% task: Pick the task:

task = 'trial_conversion';

%%%% Pick the number of:

%out: dimension of output
%in: dimension of input
Nout = 2; 
Nin = 4;  
histogram_ISI='delay';
%         histogram_ISI='target';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Tipo de Hint:         
hint='nohint'; %% No hay hint
% %         hint='hint_decision'; %% Hint de decisión

%Nhint: dimension of hint:
Nhint=1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ¿Hacemos la etapa demean?:
demean='demean'; 
% demean='nodemean';

%n: Number of s's targets to use to construct z:
% n = round(N/2);% n = round(0.4 * N); Tomar n=N sería muy costoso.
n = round(0.4 * N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters

%fracPCA: used to compute number of PCs to use
fracPCA = 0.999;

%DTRLS: number of time steps between RLS updates
%dt: timestep
%taux: decay time for rate net
%tauf: decay for fast synapses
%taus: decay for slow synapses
%taum: decay for voltage in spiking net
%Vth: spike threshold
%Vr: reset voltage

DTRLS = 2; 
dt = 0.0001;
taux = 1e-2;    %%  10 ms
taus = 1e-1;    %% 100 ms % filtro sinaptico lento
tauf = 5e-3;    %% 5 ms    % filtro sinaptico rapido 
taum = 1e-2;    %% 10 ms % membrana 
Vth = 1e-16;    %% 
Vr = -10;       %% en mV

%precompute for doing Euler integration
etaux = exp(-dt/taux);
etaus = exp(-dt/taus);
etauf = exp(-dt/tauf);
etaum = exp(-dt/taum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random connections

%Jtilde: recurrent firing rate connections
%ufout: input connections into rate network for fout
%ufin: input connecdtions into rate network for fin
%uJ: connections from low D learned signals into spiking network
%Jmu: global recurrent inhibition
%Jf: fast random synapses
%Js: slow random synapses
%s_train: which synapses will be trained in the spiking network

if any(strcmp(doing_train,{'YES'})) %String comparison; any poco útil aquí --> array de 1 elemento
    %¿Incluimos el entrenamientos? --> Sí, entonces inicializamos variables.

    randseed = 66666; %reproducir el mismo resultado a pesar de la aleatoriedad. Identidad de la red estudiada.
    rng(randseed); %set ran seed, for rand matrices

    Jtilde = gtilde * 1/sqrt(Ntilde) * randn(Ntilde);
    u_tilde_fout  = gfout * (-1 + 2 * rand(Ntilde,Nout));
    u_tilde_fin = gfin * (-1 + 2 * rand(Ntilde,Nin));
    u_tilde_fhint = gfhint * (-1 + 2 * rand(Ntilde,Nhint));

    uJ = gz * (orth((-1 + 2 * rand(N,Ntilde))));

    Jmu = muf * (1/tauf) * 1/(N) * ones(N);
    Jf = gf * (1/tauf) * 1/sqrt(N) * randn(N);
    Js = gs * (1/taus) * 1/sqrt(N) * randn(N);
    %todas las N neuronas están conectadas con todas (matrices NxN), pero
    %no todas las corrientes son no nulas.

    s_train = false(1,N);
    s_train(1,randsample(N,n)) = true;


    cd(result_archivos_folder)

        save('randseed','randseed') %save('nombre','variable')
        save('Jtilde','Jtilde')
        save('u_tilde_fout','u_tilde_fout')
        save('u_tilde_fin','u_tilde_fin')
        save('u_tilde_fhint','u_tilde_fhint')
        save('uJ','uJ')
        save('Jf','Jf')
        save('Js','Js')
        save('Jmu','Jmu')        
        save('s_train','s_train')

    cd ..

    
elseif any(strcmp(doing_train,{'NO'})) %Ya tenemos la red entrenada, sólo quiero testear.
    cd(result_archivos_folder)

        load('randseed','randseed')
        load('Jtilde','Jtilde')
        load('u_tilde_fout','u_tilde_fout')
        load('u_tilde_fin','u_tilde_fin')
        load('u_tilde_fhint','u_tilde_fhint')
        load('uJ','uJ')
        load('Jf','Jf')
        load('Js','Js')
        load('Jmu','Jmu')                
        load('s_train','s_train')    
    cd ..
end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% En la función run_model_test se generan los archivos de datos que luego
% entregamos tanto al código de phase space como al código de información mutua).
% En estos archivos de datos están ordenados los ensayos en función de las 
% clases y no en función del orden de ensayo (siguiendo así el ejemplo de los 
% datos experimentales). En otras palabras, primero recogemos aquellos ensayos en los 
% que la clase fue la número 1, a continuación, aquellos ensayos con la clase 
% número 2, y así sucesivamente. 

% 1) En el array "A" donde vamos a guardar los datos necesarios (clases, nº ensayo, 
% outcome, ...) vamos a sobre dimensionarlo. Si trabajamos con 4000 trials,
% siempre lo inicializaremos  con un número de filas mayor que el nº de
% trials. ¿Cómo lo hacemos? teniendo en cuenta tanto el número de ensayos test
% como el número de clases que tiene el set en dicha etapa test. Vamos a
% calcular cuántos ensayos se corresponderían con cada clase
% (aproximación):

    %Si el número de trials test es igual a 4000 y el número de clases es 4 (LL,SS, LS, SL).
    % 4000 trials / 4 clases = 1000 trials/clase. Entonces space_class >> 1000
        space_class=1200;

% La variable space_class nos indicará el tamaño del array "A" que le
% corresponde a una clase. Por lo que el array "A" tendrá una dimensión
% (número de filas) igual a space_class*n_clases. Como ya hemos dicho, este
% array esta sobredimensionado en un inicio (num_trials << num_filas_A).
% Esto se debe a que a medida que vayamos recorriendo los ensayos iremos
% guardando, uno a uno, los datos y la información de cada trial en este
% array en función del índice de clase. Si en el trial número 1 sale la
% clase 7, entonces nos iremos a la parte del array "A" reservado a la
% clase 7 y colocaremos ahí la información del trial. De esta manera vamos
% guardando los ensayos ordenados según su índice de clase a medida que van
% apareciendo.

% Al acabar la etapa de test sólo tendremos que quitar aquellos índices de "A"
% que no hayan sido usados, es decir, quitaremos los "huecos" extra del
% array quedándonos solo con los n_trials (4000) ya ordenados en función
% del índice de clase.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


 %% Compute PCs (Principal Components) of rate model (no se usa) 
uPC=eye(N);
% % %         uPC = run_model_PCA('PCA',TPCA,task,hint);
% Calcula componentes principales (no lo usaré)

 
%% demean spiking inputs --> Primera etapa
if any(strcmp(doing_train,{'YES'}))

    if any(strcmp(demean,{'demean'}))
        [Imu,mFR] = run_model_demean('demean',Tdemean,task,delay_values,histogram_ISI,uPC,hint);
    elseif any(strcmp(demean,{'nodemean'}))
        Imu=zeros(N,1); %Corriente de la red de espigas, para estabilizar la red
    end

elseif any(strcmp(doing_train,{'NO'}))
    cd(result_archivos_folder)
        load('Imu','Imu')
    cd ..
end


%% Train with RLS


if any(strcmp(doing_train,{'YES'}))

    [dJ_PC,w,P] = run_model_train('RLStrain',TRLS,task,delay_values,histogram_ISI,uPC,Imu,hint);

elseif any(strcmp(doing_train,{'NO'}))
    cd(result_archivos_folder)
    % Descargamos las variables resultantes del entrenamiento

        archivos_dJ_PC = dir('dJ_PC*'); % Pesos recurrentes internos de la red de espigas (ya entrenados)
        archivos_P = dir('P*'); % Matriz de correlación de las corrientes internas (ec. 35 Parga, Serrano)
        archivos_w = dir('w*'); % Pesos output de la red de espigas (ya entrenados)

        listaArchivos_dJ_PC = {archivos_dJ_PC.name}';  %% este campo contiene el nombre de los archivos
        listaArchivos_P = {archivos_P.name}';  %% este campo contiene el nombre de los archivos
        listaArchivos_w = {archivos_w.name}';  %% este campo contiene el nombre de los archivos
        clear archivos_dJ_PC; clear archivos_P; clear archivos_w;

        numArchivos = size(listaArchivos_dJ_PC,1); % Número de valores de delay que tenemos
                                                   % Tenemos tres archivos
                                                   % que descargar (dJ_PC, w, P)
                                                   % por cada delay.                                                                                                   % delay.

        for i_archivos=1:numArchivos
            load (listaArchivos_dJ_PC{i_archivos,1})
            load (listaArchivos_P{i_archivos,1})
            load (listaArchivos_w{i_archivos,1})

        end    
        clear i_archivos; clear numArchivos;
        clear listaArchivos_dJ_PC; clear listaArchivos_P; clear listaArchivos_w;
    cd ..
end

    
%for i=1:length(delay_values) % valor del delay de la tarea concreta. Si tenemos varios valores de 
    % delay, lenght trabaja con los n valores del array.
    %% Test RLS solution
    ERR_RLS = run_model_test('RLStest',Ttest,task,delay_values,histogram_ISI,uPC,Imu,dJ_PC,w,hint);

%end
