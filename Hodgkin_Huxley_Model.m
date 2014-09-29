% Aidan Dowdle
% 9/27/14
% ECE4784 Modeling Project
% Phase I: Hodgkin Huxley Model

%==========================================================================
% Function name: Hodgkin_Huxley_Model
%
% Function description: 
% This function simulates the Hodgkin Huxley model for a variety of  
% injection currents.
%==========================================================================
function Hodgkin_Huxley_Model
close all;

dt = 0.01; % Time step (ms)
% 1. Simulate a steady-state neuron with the resting membrane potential. 
%    Input current is zero.
I_inj = zeros(100/dt, 1);
run_model(I_inj, dt);

% 2. Stimulate the model neuron with a step pulse of 5 uA/cm^2 for 0.5 ms
I_inj = [5*ones(1,0.5/dt + 1) zeros(1, 99.5/dt)];
run_model(I_inj, dt);

% 3. Stimulate the model neuron with a constant current of 5 uA/cm^2
I_inj = 5*ones(1,100/dt+1);
run_model(I_inj, dt);

end

%==========================================================================
% Function name: run_model
% Inputs (2): -I_inj, the injection current array in uA/cm^2
%             -dt, the time step in ms
%
% Function description: 
% This function runs the Hodgkin Huxley model for a given injection current 
% and time step over 100 ms. The membrane voltage, K+ conductance, and Na+
% conductance are plotted afterwards.
%==========================================================================
function run_model(I_inj, dt)
% Initialize time vector and time step
t = 0:dt:100; %Time (ms)

% Initialize m, n, and h, which represent K+ channel activation, Na+
% channel activation, and Na+ channel inactivation, respectively.
m = zeros(1, length(t));
n = zeros(1, length(t));
h = zeros(1, length(t));

% Initialize alphas and betas, the rate constants
alpha_m = zeros(1,length(t));
beta_m = zeros(1,length(t));
alpha_n = zeros(1,length(t));
beta_n = zeros(1,length(t));
alpha_h = zeros(1,length(t));
beta_h = zeros(1,length(t));

% Initialize membrane voltage and conductance vectors
V_m = zeros(1,length(t));
g_K = zeros(1,length(t));
g_Na = zeros(1, length(t));

% Constants
g_K_max = 36;       %K+ channel max conductance (mS/cm^2)
g_Na_max = 120;     %Na+ channel max conductance (mS/cm^2)
g_L = 0.3;      %Leakage current conductance (mS/cm^2)
E_K = -12;      %K+ channel Nernst potential (mV)
E_Na = 115;     %Na+ channel Nernst potential (mV)
E_L = 10.6;     %Leakage channel Nernst potential (mV)
V_rest = -70;   %Resting voltage (mV)
C_m = 1.0;      %Membrance capacitance (uF/cm^2)

for i = 1:(length(t)-1)
    % Calculate all alpha and beta values
    alpha_m(i) = 0.1*((25 - V_m(i))/(exp((25-V_m(i))/10)-1));
    beta_m(i) = 4*exp(-1*V_m(i)/18);
    alpha_n(i) = 0.01*((10-V_m(i))/(exp((10-V_m(i))/10)-1));
    beta_n(i) = 0.125*exp(-1*V_m(i)/80);
    alpha_h(i) = 0.07*exp(-1*V_m(i)/20);
    beta_h(i) = 1/(exp((30-V_m(i))/10)+1);
    
    % Initial conditions
    if(i == 1)
        m(i) = alpha_m(i)/(alpha_m(i) + beta_m(i));
        n(i) = alpha_n(i)/(alpha_n(i) + beta_n(i));
        h(i) = alpha_h(i)/(alpha_h(i) + beta_h(i));
    end
    
    % Calculate conductances
    g_Na(i) = m(i)^3*g_Na_max*h(i);
    g_K(i) = n(i)^4*g_K_max;
    
    % Calculate currents
    I_Na = g_Na(i)*(V_m(i) - E_Na);
    I_K = g_K(i)*(V_m(i) - E_K);
    I_L = g_L*(V_m(i) - E_L);
    I_ion = I_inj(i) - I_K - I_Na - I_L;
    
    % Calculate membrane voltage, m, n, and h using Euler's method
    V_m(i+1) = V_m(i) + I_ion/C_m*dt;
    m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*dt;
    n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*dt;
    h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*dt;
end

V_m = V_m + V_rest;

% Plot membrane potential
figure;
plot(t, V_m)
axis([0, 100, -100, 40])
title('Membrane Potential')
xlabel('Time (ms)')
ylabel('Voltage (mV)')

% Plot conductances
figure;
plot(t, g_Na, 'r', t, g_K, 'b');
title('gK and gNa')
xlabel('Time (ms)')
ylabel('Conductance (S/cm2)')
end