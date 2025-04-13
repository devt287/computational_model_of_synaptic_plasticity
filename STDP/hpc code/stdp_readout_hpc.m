task_id = str2double(getenv('SLURM_ARRAY_TASK_ID'));
a=0.25;
b=0.9;
pre=1.3;
% Parameters
K0=0.5; k1=2/10; k2=15/10; k3=1/10; k4=a*120/10; k11=2/10; k12=15/10; k13=1/10; 
k14=a*80/10;
Km1=10; Km2=0.3;  Km12=b; 
Km=3;
Km4=4;
Ktot=20;Ptot=20;P0=0.5;Ca_basal=0.1;
c_1=1; c_2=1; c_3=6; c_4=8; tau_glut_f=50;
tau_glut_s=200; I_glut_f=0.5; I_glut_s=0.5; V_Ca=130; MG2=1.5; V_VDCC=130; I_bpap_f=0.75; I_bpap_s=0.25;
tau_bpap_f=3; tau_bpap_s=25; Vmax=100; tau_decay=12;
timespan=150000;tau_1_ep=50; tau_2_ep=5; s=1; span = [0 timespan];
v=0.5;
Km11=10;



%a is the parameter that apply to k4 and k14 lowering down the sensitivity
%to the basal state.
ratio=100;
pK_init=0.018;
P_init=0.085;
C_init=5.708;
C_tot=(k11*(Ptot-P_init)*P_init)/(Km11+Ptot-P_init);  
Km11=Ptot/ratio;

%k11_new=C_tot*(Km11+Ptot-P_init)/((Ptot-P_init)*P_init);
k11_new=0.045;
init_val=[pK_init,P_init,C_init,Ca_basal];


% Define delta_t values
%delta_t_values = [-500, -400, -300, -200, -150, -100, -75, -50, -25, -10, 0, 10, 25, 50, 75, 100, 200, 300, 400, 500];
delta_t_values = [-500, -400, -300, -200, -150, -100, -75, -50, -25, -20,-15,-10, 0, 10,12,14,16,18,20,22 25, 50, 75, 100, 200, 300, 400, 500];
fid = fopen('Km=3_K14=2_k11=0.045_1.3_fine.txt', 'a');

G_VDCCs = 0.01:0.01:2.5;
G_NMDAs=G_VDCCs*1.5;


G_VDCC = G_VDCCs(task_id)/50;  % Use only the specific G_VDCC for this job
G_NMDA = G_NMDAs(task_id)/50 ;  % Corresponding G_NMDA

result_array = zeros(1, length(delta_t_values));
for idx = 1:length(delta_t_values)    
    %par2=[beta,lambda]
    par2=[0.9,1/5];
    %[index,thr,h]
    mu_par = [1,2,40];
    % Calculate t_bpap
    t_glut=1050;
    t_bpap = t_glut + delta_t_values(idx);
    t_glut_list = t_glut + 1000 * (0:59);
    t_bpap_list = t_bpap + 1000 * (0:59);
    
    theta1 = @(t) sum((t>t_glut_list(t_glut_list<t)).*(t<t_glut_list(t_glut_list<t)+300));
    theta2 = @(V,t) sum((V>-30).*(t>t_bpap_list(t_bpap_list<t)).*(t<t_bpap_list(t_bpap_list<t)+2));
    H=@(V) (V_Ca-V)/(1+(MG2/3.57)*exp(-0.062*V)); 
    I_NMDA = @(V,t) sum(P0*G_NMDA.*theta1(t).*(I_glut_f*exp((t_glut_list(t_glut_list<t)-t)/tau_glut_f)+I_glut_s*exp((t_glut_list(t_glut_list<t)-t)/tau_glut_s)).*H(V));
    I_VDCC= @(V,t) G_VDCC.*theta2(V,t).*(V_VDCC-V);
    V_bpap = @(t) Vmax*sum(((I_bpap_f*exp(-(t-t_bpap_list(t_bpap_list<t))/tau_bpap_f))+(I_bpap_s*exp(-(t-t_bpap_list(t_bpap_list<t))/tau_bpap_s))));
    % Compute norm
    t_values = linspace(min(span), max(span), timespan);
    fun = exp(-(t_values - t_glut) / tau_1_ep) - exp(-(t_values - t_glut) / tau_2_ep);
    norm = max(fun);
    % Now that norm is defined, define V_Epsp
    V_Epsp = @(t) (s/norm)*sum((exp(-(t-t_glut_list(t_glut_list<t))/tau_1_ep)-exp(-(t-t_glut_list(t_glut_list<t))/tau_2_ep)));
    % Compute V values
    V = @(t) V_bpap(t) + V_Epsp(t) - 65;
    g=@(t,Ca) I_NMDA(V(t), t)+I_VDCC(V(t), t)-(Ca-Ca_basal)/tau_decay; 
    % Solve the differential equation        
    opts = odeset('RelTol',1e-7, 'AbsTol',1e-6, 'MaxStep',0.1); 
    mu = @(Ca) (Ca^mu_par(1))./((Ca^mu_par(1)) + mu_par(2)^mu_par(1))*mu_par(3);
    nu = @(pK) v*par2(1)/(1+ par2(2)*pK)+1-par2(1);
    % Differential equations for Y(1) and Y(2)
    g2 = @(t,Y,Ca)[k1*((Ktot-Y(1)-Y(3))/(Km1+(Ktot-Y(1)-Y(3))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca.^4)*(Ktot-Y(1)-Y(3)))/(Km4^4+Ca.^4);
    (k11_new*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(pre*Y(1)+K0)+k13*P0+(k14*(Ca.^3)/(Km^3+Ca.^3))*(Ptot-Y(2)));
    -mu(Ca)*Y(3)+nu(Y(1))*(Ktot-Y(1)-Y(3));];
    g_tot = @(t,Y) [g2(t,Y(1:3), Y(4))/1000; g(t,Y(4))]; 
    % Solve the differential equations for Y(1) and Y(2)
    [t1,X]=ode45(g_tot, span, init_val, opts);
    difference = X(end,1) - X(end,2);
    if difference > 1
        result_array(idx) = 1;
    elseif difference < -1
        result_array(idx) = -1;
    else
        result_array(idx) = 0;
    end
end
fprintf(fid, 'G_VDCC=%.2f,G_NMDA=%.2f, [%s]\n', G_VDCC*50, G_NMDA*50, strjoin(arrayfun(@num2str, result_array, 'UniformOutput', false), ','));

fclose(fid);