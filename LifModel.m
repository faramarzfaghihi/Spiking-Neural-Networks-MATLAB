%%%%%MATLAB code for adaptive LIF neuron written by Faramarz Faghihi%%%%

clear

% membrane constants
tau = 0.030;
R = 3e7;
E = 0;
thresh = 0.035;

% Time
dt = 0.001;
T = 0.5;
No_steps = T ./ dt;
ts = linspace(0, T, No_steps + 1);

% Injection current
I_0 = 2e-9;
I_len = 0.25;
I_start = 0.05;
I_start_index = I_start ./ dt
I_finish_index = (I_start + I_len) ./ dt;
I = zeros(1, No_steps + 1);
I(I_start_index:I_finish_index) = I_0;

% membrane potential
V_0 = 0;
V = zeros(1, No_steps + 1);
V(1) = V_0;
t_spikes = [];
%-------------------------------------------------------------------------
for i=1:No_steps
    dV = (1 ./ tau) .* (E - V(i) + I(i) .* R) .* dt;
    V(i+1) = V(i) + dV;
    if V(i+1) > thresh
        V(i+1) = E;
        t_spikes = [t_spikes (i - 1) * dt];
    end
end
%-------------------------------------------------------------------------

spike_height = 0.1;
No_spikes = length(t_spikes);

if No_spikes > 0
    spts = [t_spikes; t_spikes];
    y1 = thresh .* ones(1, No_spikes);
    y2 = y1 + spike_height;
    sp = [y2; y1];
end

%--------------------------------------------------------------------------

clf
plot(ts, V);
xlabel('Time') 
ylabel('Membrane potential') 
title('Membrane potential of single neurons')
hold on
if No_spikes > 0
    h = line(spts, sp);
    for k=1:length(h)
        set(h(k), 'Color', [0 0 1])
    end
end

hold off
