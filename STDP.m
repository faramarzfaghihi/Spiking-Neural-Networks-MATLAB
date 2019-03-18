%%%%%%% MATLAB code for STDP rule in Diehl et al., 2015 
%%%%   written by Faramarz Faghihi

close all;

fileID = fopen('t10k-images.idx3-ubyte', 'r');

A = fread(fileID, 1, 'uint32');
magicNumber = swapbytes(uint32(A));

A = fread(fileID, 1, 'uint32');
totalImages = swapbytes(uint32(A));

A = fread(fileID, 1, 'uint32');
numRows = swapbytes(uint32(A));

A = fread(fileID, 1, 'uint32');
numCols = swapbytes(uint32(A));

for k = 1 : totalImages
    
    A = fread(fileID, numRows*numCols, 'uint8');
    imageCellArray{k} = reshape(uint8(A), numCols, numRows)';
    
end

%//Close the file
fclose(fileID);

image1 = imageCellArray{1, 80};

dbinput = double(image1)/400;
input_vec = reshape(dbinput.',1,[]); % image as 784 lenght vector
Input_spike_trains = [];

for t =1:length(input_vec)
    
    Input_spike_trains(t,:) = SpikeGen(input_vec(t));
    
end

%// spike_trains : spike trains of each 784 input neurons

%----------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%% parameters value of IF model %%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(Input_spike_trains);   % size of input
Exc_spike_hist = zeros (1,N); % initial spike history of one excitatory neurons
V_th = -52;                   % mv   threshold of firing of excitatoy neurons
V_reset = -65;                % Reset value
EL = -65.3;                   % [mV] leak reversal potential
E_exc = -65;                  % Equilibrium potentials of excitatorY synapses
Vm = zeros (1,N);             % Initial value of membrane potential of excitatory neuron
Vm(1) = -25;

%%%%%%%%%%%%%%%%%%% parameters of synapse dynamics %%%%%%%%%%%%%%%%%%%%%%

X_pre = zeros(M,N);           % Initial values of presynaptic traces
dt = 0.001;                   % Time step used in IF model
ta_g = 1.2;                   % Decay rate of excitatory conductance  
ta = 0.005;                   % Decay rate of membrane potential
ta_trace = 40;                % Decay rate of presynaptic trace
Weight = zeros (M,N)+0.05;           % Intial weight matrix
weightFinal = zeros (1,N);           % Initial total induced weight change between 
                                     %--> one excitatory neuron and input layer
G_e = zeros (M,N);                   % Initial values of excitatory conductance

%%%%%%%%%%%%%% synaptic conductance change induced by input layer %%%%%%%%%

for  i = 1:M    % loop for neurons in input lyer = 784
    for j = 2:N % loop for time
        if Input_spike_trains(i,j)>0
            G_e(i,j) = G_e(i,j-1)- G_e(i,j-1)/ta_g + 1;
        else
            G_e(i,j) = G_e(i,j-1)- G_e(i,j-1)/ta_g;
        end
    end
end

G_e_Exc = sum(G_e,1); % Sum of change in g_e induced by presynaptic spikes


%%%%%%%%%%%%%%%%% presynaptic trace dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for  i = 1:M         % loop for neurons in input lyer = 784
    for j = 2:N      % loop for time
        
        if Input_spike_trains(i,j)>0
            X_pre(i,j) = X_pre(i,j-1)- X_pre(i,j-1)/ta_trace + 1;
        else
            X_pre(i,j) = X_pre(i,j-1)- X_pre(i,j-1)/ta_trace;
            
        end
    end
    
end

%%%%%%%%%%%%%%% STDP learning parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ro =0.05;    % dependenc parametere of the update on the previous weight
mu = 1.8;    % Learning rate
X_tar = 2;   % target value of presynaptic trace at the moment of a post synaptic spike
W_max = 2.5; % Maximum weight

%------------------------------------------------------------------------


for t=2:N                                     % loop for time
    
    spike = (Input_spike_trains(:,t));        % Input into excitatory neuron at each time bin
    
    %%%%%%%%%%%%%%%%% Integrate and fire model of excitatory neuron %%%%%%%%%%%%%%%%%
    
    Vm(t) = Vm(t-1) +(1/ta)*dt*( (EL -Vm(t-1)) + G_e_Exc (1,t)*(E_exc -Vm(t-1))); % Membrane potential equation
    
    if Vm(t) > V_th
        Exc_spike_hist (t)=1;
        Vm(t) = V_reset;
    end
    %%%%%%%%%%%% STDP Learning %%%%%%%%%%%%%%%%%%
    if Exc_spike_hist (t)==1
        
        for k = 1:M                            % loop for neurons in input lyer = 784
              
            Weight(k,t) =  Weight(k,t-1)+ ro* (X_pre(k,t)-X_tar)*((W_max-Weight(k,t-1))) ^mu;
            if Weight(k,t)<0
                Weight(k,t)=0;
            end
        end
    end
    
    if Exc_spike_hist (t)==0
       % 'no change in synaptic weights' 
        for k = 1:M                             % loop for neurons in input lyer = 784
            
            Weight(k,t) =  Weight(k,t-1);
            if Weight(k,t)<0
                Weight(k,t)=0;
            end
        end
    end
    weightFinal(1,t) = sum(Weight(:,t));
    
end

%--------------------------------------------------------------------------
figure
plot(Exc_spike_hist)
title('Excitatory neuron spiking and spike train of a excitatory neuron')
xlabel ('timebins')
ylabel ('pre-synaptic trace')
hold on
plot(Input_spike_trains(521,:)) % spike history of a neuron in input layer
hold on
plot(X_pre(521,:),'g')    % pre synaptic trace of a neuron in input layer
legend({'ExcSpike','inputSpike'},'Location','northwest')

figure
plot(Weight(521,:),'r')   % synaptic weight of a neuron in input layer
title('Synaptic weight chnage of a excitatory neuron')
xlabel ('timebins')
ylabel ('synaptic weight')
fir_Rate_Exc_n = mean (Input_spike_trains(521,:))  % mean of firing rate of excitatory neurons



