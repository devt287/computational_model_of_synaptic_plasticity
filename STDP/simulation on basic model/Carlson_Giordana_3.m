% Parameters
K0=0.5; k1=2; k2=15; k3=1; k4=120; k11=2; k12=15; k13=1; k14=80; %1/s
Km1=10; Km2=0.3; Km11=10; Km12=1; Km=4; Ktot=20; Ptot=20; Atot=1;
c_1=1; c_2=1; c_3=6; c_4=8; P0=0.5; G_NMDA=1.75/50; tau_glut_f=50;
tau_glut_s=200; I_glut_f=0.5; I_glut_s=0.5; V_Ca=130; MG2=1.5;
G_VDCC=0.8/50; t_glut=1050; V_VDCC=130; I_bpap_f=0.75; I_bpap_s=0.25;
tau_bpap_f=3; tau_bpap_s=25; Vmax=100; tau_decay=12; Ca_basal=0.1;
tau_1_ep=50; tau_2_ep=5; s=1; span = [0 1500];

% Initial Ca value
Ca_init = 0.1;

% Define delta_t values
delta_t_values = [-45, -30, -15, 0, 15, 30, 45];

% Figure initialization
figure;
hold on;

for idx = 1:length(delta_t_values)
    % Calculate t_bpap
    t_bpap = t_glut + delta_t_values(idx);
    % Define auxiliary functions
    theta1 = @(t) (t>t_glut);
    theta2 = @(V,t) (V>-30).*(t>t_bpap).*(t<t_bpap+2);
    H=@(V) (V_Ca-V)/(1+(MG2/3.57)*exp(-0.062*V)); 
    I_NMDA = @(V,t) P0*G_NMDA.*theta1(t).*(I_glut_f*exp((t_glut-t)/tau_glut_f)+I_glut_s*exp((t_glut-t)/tau_glut_s)).*H(V);
    I_VDCC= @(V,t) G_VDCC.*theta2(V,t).*(V_VDCC-V);
    V_bpap = @(t) Vmax*((I_bpap_f*exp(-(t-t_bpap)/tau_bpap_f))+(I_bpap_s*exp(-(t-t_bpap)/tau_bpap_s)))*(t>t_bpap);

    % Compute norm
    t_values = linspace(min(span), max(span), 10000);
    fun = exp(-(t_values - t_glut) / tau_1_ep) - exp(-(t_values - t_glut) / tau_2_ep);
    norm = max(fun);

    % Now that norm is defined, define V_Epsp
    V_Epsp = @(t) (s/norm)*(exp(-(t-t_glut)/tau_1_ep)-exp(-(t-t_glut)/tau_2_ep))*(t>t_glut);

    % Compute V values

    V = @(t) V_bpap(t) + V_Epsp(t) - 65; %NEW

    % Define the main function for ode45
    g=@(t,Ca) I_NMDA(V(t), t)+I_VDCC(V(t), t)-(Ca-Ca_basal)/tau_decay; %NEW

    % Solve the differential equation
    opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',0.1); 

    % Differential equations for Y(1) and Y(2)
    g2 = @(t,Y,Ca)[k1*((Ktot-Y(1))/(Km1+(Ktot-Y(1))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca.^4)*(Ktot-Y(1)))/(Km^4+Ca.^4);
    (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca.^3)/(Km^3+Ca.^3))*(Ptot-Y(2)))];
    
    g_tot = @(t,Y) [g2(t,Y(1:2), Y(3))/1000; g(t,Y(3))]; 

    % Solve the differential equations for Y(1) and Y(2)
    [t1,X]=ode45(g_tot, span, [0;0;0.1], opts); 

    % Plot results
    
    plot(t1, X(:,1), 'LineWidth', 1.5, 'DisplayName', ['pK, \Delta_t = ', num2str(delta_t_values(idx))]);
    xlim([1000 1200]);
end
xlabel('Time');
ylabel('pK');
legend('Location', 'best');
hold off;