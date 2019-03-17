%%%%%% MATLAB code for excitatory inhibitory neural network using Izhikevich model%%%%%%%%
%%%%%%

close all
clear all
clc

number_INH_Neurons = 95;
number_EXC_Neurons = 95;
number_OUTPUT_Neurons = 30;
time_stimulation = 4000;
m1 = number_INH_Neurons;
m2 = number_EXC_Neurons;
m3 = number_OUTPUT_Neurons;
nt = time_stimulation;
SP_INH = zeros(m1,nt);
SP_EXC = zeros(m2,nt);
firing_i = 0.2;
firing_e = 0.8;
firing_prob_INH = firing_i;
weigth_EXC = 0.5;
connectivityRate_INH = 0.9;
connectivityRate_EXC = 0.9;
initial_weigth = 6;
firing_prob_EXC = firing_e;
ta23 = 3;
ta_inh = 12;

%%%%%%%%  Inhibitory networks    %%%%%%%%

for k1 = 1:m1
    
    firing_prob_INH = firing_i;
    
    ran1 = rand (1,nt);
    
end

SP_INH(k1,:) = ran1 <= firing_prob_INH;


%%%%%%%%%%%%%%% Connectivity%%%%%%%%%%%%%%%%%%

for k1 = 1:m1
    
    for k2 = 1:m3
        
        ran1 = rand (1,m3);
        
    end
    
    ConnINH_OUT(k1,:) = ran1 <= connectivityRate_INH;
    
end


%%%%%%%% Excitatory networks     %%%%%%%%


for k1 = 1:m2
    
    ran1 = rand (1,nt);
    SP_EXC(k1,:) = ran1 <= firing_prob_EXC;
    
end

%%%%%%%%%%%%%%% Connectivity matric %%%%%%%%%%%%%%%%%%

for k1 = 1:m2
    
    for k2 = 1:m3
        
        ran1 = rand (1,m3);
        
    end
    
    ConnEXC_OUT(k1,:) = ran1 <= connectivityRate_EXC;
    
end

spike_EXC = weigth_EXC*SP_EXC'*ConnEXC_OUT;

input_spike_EXC = spike_EXC';

spike_INH = SP_INH'*ConnINH_OUT;
input_spike_INH = spike_INH';

weigth_INH = initial_weigth*ones(m3,nt); %%%%%% Inhibitory neurons weigth matrics %%%%%%

%%%%%%%% Sparse spiking in DG network  %%%%%%%%%%

C =100;
vr =-60;
vt =-40;
k =0.7;
a =0.03; b =-2;c =-50;d =100;
vpeak =35;
T =nt; tau =0.9;
n =round(T/tau);
v =vr*ones(m3,n); u=0;
dt = 0.01;

%--------------------------------------------------------------------------

for t = 1:m3
    
    activator = 0;
    
    for i = 1:n-1
        
        if i<nt
            %%%%%%%%%%%% updates of weigths m2,m3 %%%%%%%%%%%%%%%%
            
            for tt = 1:100
                
                RM (tt) = activator*(tt/ta23) * exp(-tt/ta23);
                
            end
            
            stim = sum(RM);
            
            weigth_INH (t,i+1) = weigth_INH (t,i)-dt*(1/ta_inh)*(weigth_INH (t,i))+stim;
            
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
            I = input_spike_EXC(t,i)-weigth_INH(t,i+1)*input_spike_INH(t,i);
            v(t,i+1 ) = v(t,i)+tau*(k*(v(t,i)-vr)*(v(t,i)-vr)-u(i)+I)/C;
            u(i+1) = u(i)+tau*a*(b*(v(t,i)-vr)-u(i));
            
            if v(t,i+1)>= vpeak
                activator = 1;
            else
                activator = 0;
            end
            if v(t,i+1)>= vpeak
                v(t,i) = vpeak;
                v(t,i+1) = c;
                u(i+1) = u(i+1)+d;
                
            end
            
        end
        
    end
    
end

%---------------------------------------------------------------------------
figure
plot(tau*(1:n),v(1,:));
xlabel('Time') 
ylabel('Membrane potential') 
title('Membrane potential changes of a single neuron')
