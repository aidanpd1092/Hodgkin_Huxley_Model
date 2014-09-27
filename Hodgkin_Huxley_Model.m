function Hodgkin_Huxley_Model
% Time vector
dt = 0.01; %ms
t = 0:dt:100 %ms

m = zeros(1, length(t));
n = zeros(1, length(t));
h = zeros(1, length(t));

alpha_m = zeros(1,length(t));
beta_m = zeros(1,length(t));
alpha_n = zeros(1,length(t));
beta_n = zeros(1,length(t));
alpha_h = zeros(1,length(t));
beta_h = zeros(1,length(t));

V_m = zeros(1,length(t));

% Constants
g_K = 36; %mS/cm^2
g_Na = 120; %mS/cm^2
g_L = 0.3; %mS/cm^2
E_K = -12; %mV
E_Na = 115; %mV
E_L = 10.6; %mV
V_rest = -70; %mV
C_m = 1.0; %uF/cm^2
I = 0;

V_m(1) = -70; %mV
alpha_m(1) = 0.1*((25 - V_m(1))/(exp((25-V_m(1))/10)-1));
beta_m(1) = 4*exp(-1*V_m(1)/18);
alpha_n(1) = 0.01*((10-V_m(1))/(exp((10-V_m(1))/10)-1));
beta_n(1) = 0.125*exp(-1*V_m(1)/80);
alpha_h(1) = 0.07*exp(-1*V_m(1)/20);
beta_h(1) = 1/(exp((30-V_m(1))/10)+1);
m(1) = alpha_m(1)/(alpha_m(1) + beta_m(1));
n(1) = alpha_n(1)/(alpha_n(1) + beta_n(1));
h(1) = alpha_h(1)/(alpha_h(1) + beta_h(1));

for i = 1:(length(t)-1)
    I_Na = m(i)^3*g_Na*h(i)*(V_m(i) - E_Na);
    I_K = n(i)^4*g_K*(V_m(i) - E_K);
    I_L = g_L*(V_m(i) - E_L);
    I_ion = I - I_K - I_Na - I_L;
    
    V_m(i+1) = V_m(i) + I_ion/C_m*dt;
    m(i+1) = m(i) + (alpha_m(i)*(1-m(i))+beta_m(i)*m(i))*dt;
    n(i+1) = n(i) + (alpha_n(i)*(1-n(i))-beta_n(i)*n(i))*dt;
    h(i+1) = h(i) + (alpha_h(i)*(1-h(i))-beta_h(i)*h(i))*dt;
    
    alpha_m(i+1) = 0.1*((25 - V_m(i+1))/(exp((25-V_m(i+1))/10)-1));
    beta_m(i+1) = 4*exp(-1*V_m(i+1)/18);
    alpha_n(i+1) = 0.01*((10-V_m(i+1))/(exp((10-V_m(i+1))/10)-1));
    beta_n(i+1) = 0.125*exp(-1*V_m(i+1)/80);
    alpha_h(i+1) = 0.07*exp(-1*V_m(i+1)/20);
    beta_h(i+1) = 1/(exp((30-V_m(i+1))/10)+1);
end

plot(t, V_m);
end