% Define the parameter value that we want to tune
parameter_pairs = [
    0.9, 1.5, 1.47, 100, 0.25, 0.0375,60,1,1.3,3;
    0.9, 1.5, 1.47, 100, 0.25, 0.040,60,1,1.3,3;
];






%Define the delta_t values
%delta_t_values=[50,75,100];
%delta_t_values = [-25,100,110,120,130,140,150,160,170,180,190 200];
%delta_t_values = [-500,-200,-100,-75, -50, -25,-10,0,10,12,16,20,22,25,50,75,100,200,500];
delta_t_values = [-500, -400, -300, -200, -150, -100, -75, -50, -25, -10, 0, 10, 25, 50, 75, 100, 200, 300, 400, 500];







parfor iter = 1:size(parameter_pairs, 1)
    G_VDCC = parameter_pairs(iter, 3) / 50;

    %pre is the prefactor defined in the equation. see the equation in the paper
    pre = parameter_pairs(iter, 9);

    a = parameter_pairs(iter, 5);  
    b = parameter_pairs(iter, 1);
    c = parameter_pairs(iter, 2);
    pairs= parameter_pairs(iter, 7)-1;
    ratio = parameter_pairs(iter, 4);
    k11 = 0.2;
    hertz=parameter_pairs(iter, 8),
    frequency_time=1000/hertz;
    
    % Define constant parameters
    K0 = 0.5; k1 = 2/10; k2 = 15/10; k3 = 1/10; k4 = a * 120 / 10; 
    k12 = 15/10; k13 = 1/10; k14 = a * 80 / 10; Km1 = 10; Km2 = 0.3; 
    Km12 = b; 
    Km=parameter_pairs(iter, 10);
    Km4=4
    Ktot = 20; Ptot = 20; P0 = 0.5; Ca_basal = 0.1;

    % Initial values
    pK_init=0.018;
    P_init=0.085;
    C_init=5.708;


    % Calculate Km11
    Km11=10;
    C_tot=(k11*(Ptot-P_init)*P_init)/(Km11+Ptot-P_init);  
    Km11=Ptot/ratio;
    
    %k11_new=C_tot*(Km11+Ptot-P_init)/((Ptot-P_init)*P_init); another definition of k11
    k11_new=parameter_pairs(iter, 6);
    
    % parameter largely followed by Carlson & Giordana's paper
    tau_glut_f = 50; tau_glut_s = 200; I_glut_f = 0.5; I_glut_s = 0.5; 
    V_Ca = 130; MG2 = 1.5; V_VDCC = 130; I_bpap_f = 0.75; I_bpap_s = 0.25; 
    tau_bpap_f = 3; tau_bpap_s = 25; Vmax = 100; tau_decay = 12; 
    timespan = 150000; tau_1_ep = 50; tau_2_ep = 5; s = 1; span = [0 timespan]; 
    v = 0.5; init_val = [0.018, 0.085, 0.2, Ca_basal];

    % Create a figure with two rows: one for pK/P/C, and one for calcium
    figure('visible', 'off');
    set(gcf, 'Position', [100, 100, 2400, 1600]); % Adjust figure height for two rows
    
    G_NMDA = G_VDCC * c;
    
    % Initialize result array
    result_array = zeros(1, length(delta_t_values));
    
    for idx = 1:length(delta_t_values)
        % Recalculate everything for each delta_t and k11
        t_glut = 1050;
        t_bpap = t_glut + delta_t_values(idx);
        t_glut_list = t_glut + frequency_time * (0:pairs);
        t_bpap_list = t_bpap + frequency_time * (0:pairs);

        theta1 = @(t) sum((t > t_glut_list(t_glut_list < t)).*(t < t_glut_list(t_glut_list < t) + 300));
        theta2 = @(V, t) sum((V > -30).*(t > t_bpap_list(t_bpap_list < t)).*(t < t_bpap_list(t_bpap_list < t) + 2));
        H = @(V) (V_Ca - V) / (1 + (MG2 / 3.57) * exp(-0.062 * V));
        I_NMDA = @(V, t) sum(P0 * G_NMDA .* theta1(t) .* (I_glut_f * exp((t_glut_list(t_glut_list < t) - t) / tau_glut_f) + I_glut_s * exp((t_glut_list(t_glut_list < t) - t) / tau_glut_s)) .* H(V));
        I_VDCC = @(V, t) G_VDCC .* theta2(V, t) .* (V_VDCC - V);
        V_bpap = @(t) Vmax * sum(((I_bpap_f * exp(-(t - t_bpap_list(t_bpap_list < t)) / tau_bpap_f)) + (I_bpap_s * exp(-(t - t_bpap_list(t_bpap_list < t)) / tau_bpap_s))));
        
        t_values = linspace(min(span), max(span), timespan);
        fun = exp(-(t_values - t_glut) / tau_1_ep) - exp(-(t_values - t_glut) / tau_2_ep);
        norm = max(fun);
        V_Epsp = @(t) (s / norm) * sum((exp(-(t - t_glut_list(t_glut_list < t)) / tau_1_ep) - exp(-(t - t_glut_list(t_glut_list < t)) / tau_2_ep)));
        V = @(t) V_bpap(t) + V_Epsp(t) - 65;
        g = @(t, Ca) I_NMDA(V(t), t) + I_VDCC(V(t), t) - (Ca - Ca_basal) / tau_decay;
        opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-6, 'MaxStep', 0.1);
        mu = @(Ca) (Ca^1)./((Ca^1) + 2^1) * 40;
        nu = @(pK) v * 0.9 / (1 + (1/5) * pK) + 1 - 0.9;

        g2 = @(t,Y,Ca)[k1*((Ktot-Y(1)-Y(3))/(Km1+(Ktot-Y(1)-Y(3))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca.^4)*(Ktot-Y(1)-Y(3)))/(Km4^4+Ca.^4);
        (k11_new*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(pre*Y(1)+K0)+k13*P0+(k14*(Ca.^3)/(Km^3+Ca.^3))*(Ptot-Y(2)));
        -mu(Ca)*Y(3)+nu(Y(1))*(Ktot-Y(1)-Y(3));];

        g_tot = @(t, Y) [g2(t, Y(1:3), Y(4)) / 1000; g(t, Y(4))];
        [t1, X] = ode45(g_tot, span, init_val, opts);
        
        % First row subplots for pK, P, C
        subplot(2, length(delta_t_values), idx); % First row for pK, P, C
        plot(t1, X(:, 1), 'LineWidth', 1.5, 'DisplayName', 'K');
        hold on;
        plot(t1, X(:, 2), 'LineWidth', 1.5, 'DisplayName', 'P');
        plot(t1, X(:, 3), 'LineWidth', 1.5, 'DisplayName', 'C');
        xlim([0 timespan]);
        set(gca, 'YScale', 'log');
        ylim([1e-3, 24]);
        legend(['\Deltat = ', num2str(delta_t_values(idx)), ' ms'], 'Location', 'best');
        
        % Second row subplots for calcium
        subplot(2, length(delta_t_values), idx + length(delta_t_values)); % Second row for calcium
        plot(t1, X(:, 4), 'LineWidth', 1.5, 'DisplayName', 'Calcium');
        xlim([0 timespan]);
        ylim([0 50]);
        legend(['Ca (\Deltat = ', num2str(delta_t_values(idx)), ' ms'], 'Location', 'best');
    end
    
    % Update the title
    sgtitle(['1:', num2str(c), ' short_G.VDCC = ', num2str(G_VDCC * 50), ', G.NMDA = ', num2str(G_NMDA * 50), ', Km12= ', num2str(Km12), ', a= ', num2str(a), ', ratio= ', num2str(int32(ratio)), ', k11= ', num2str(k11_new),", pre factor= ",num2str(pre),", Km= ",num2str(Km)]);
    filename = sprintf('1:%2f_short_G_VDCC=%.3f_G_NMDA=%.3f_Km12=%.3f_a=%.2f_ratio=%d_k11=%.3f_pre=%.2f_hertz=%.2f_Km=%.1f.png', c, G_VDCC * 50, G_NMDA * 50, Km12, a, int32(ratio), k11_new,pre,hertz,Km);
    print(filename, '-dpng');
end