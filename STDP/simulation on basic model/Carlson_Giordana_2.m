% Define parameters
P0=0.5;  G_NMDA=1.75/50; tau_glut_f=50; tau_glut_s=200; I_glut_f=0.5; I_glut_s=0.5; 
V_Ca=130; MG2=1.5; G_VDCC=0.8/50; t_glut=1050; V_VDCC=130; I_bpap_f=0.75; 
I_bpap_s=0.25; tau_bpap_f=3; tau_bpap_s=25; Vmax=100; tau_decay=12; Ca_basal=0.1;
tau_1_ep=50; tau_2_ep=5; s=1; span = [0 1500]; 

% Define delta_t values
delta_t_values = [-45, -30, -15, 0, 15, 30, 45];

figure;
hold on;

% Prepare for the legend
legend_entries = cell(length(delta_t_values), 1);

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
    [t,Ca] = ode45(g,span,0.1,opts);

    % Plot the result
    plot(t, Ca, 'LineWidth', 1.5);
    xlim([1000 1200]);

    % Add to legend entries
    legend_entries{idx} = ['delta_t = ', num2str(delta_t_values(idx))];
end
% Add labels to the axes
xlabel('Time');
ylabel('Ca');
legend(legend_entries, 'Location', 'best');
hold off;