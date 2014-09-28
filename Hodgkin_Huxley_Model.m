% Aidan Dowdle
% 9/27/14
% ECE4784 Modeling Project
% Phase I: Hodgkin Huxley Model
%
% Function Hodgkin_Huxley_Model provides a computational model for the
% Hodgkin Huxley neuron model.
function Hodgkin_Huxley_Model
dt = 0.01; % Time step (ms)
% 1. Simulate a steady-state neuron with the resting membrane potential. 
%    Input current is zero.
I = zeros(100/dt, 1);
[V_m, g_K, g_Na, t] = run_model(I, dt);
plot(t, V_m);

figure(2);
plot(t, g_K, 'b', t, g_Na, 'r');

% 2. Stimulate the model neuron with a step pulse of 5 uA/cm^2 for 0.5 ms
I = [5*ones(1,0.5/dt) zeros(1, 99.5/dt + 1)];
[V_m, g_K, g_Na, t] = run_model(I, dt);
figure(3);
plot(t, V_m);

figure(4);
plot(t, g_K, 'b', t, g_Na, 'r');

% 3. Stimulate the model neuron with a constant current of 5 uA/cm^2
I = 5*ones(1,100/dt+1);
[V_m, g_K, g_Na, t] = run_model(I, dt);
figure(5);
plot(t, V_m);

figure(6);
plot(t, g_K, 'b', t, g_Na, 'r');

end

function [V_m, g_K, g_Na, t] = run_model(I, dt)
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
    I_ion = I(i) - I_K - I_Na - I_L;
    
    % Calculate membrane voltage, m, n, and h using Euler's method
    V_m(i+1) = V_m(i) + I_ion/C_m*dt;
    m(i+1) = m(i) + (alpha_m(i)*(1-m(i))-beta_m(i)*m(i))*dt;
    n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*dt;
    h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*dt;
end

V_m = V_m + V_rest;
end