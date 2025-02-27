P0=0.5;  G_NMDA=1.75;     tau_glut_f=50;     tau_glut_s=200;    I_glut_f=0.5;      I_glut_s=0.5;    V_Ca=130;      MG2=1.5;
G_VDCC=0.8;    t_glut=1050;  t_bpap=1065;   V_VDCC=130;   I_bpap_f=0.75;    I_bpap_s=0.25;  tau_bpap_f=3;   tau_bpap_s=25;
Vmax=100;   tau_decay=12;   Ca_basal=0.1;   V=20;

span = [0 2000];

theta1 = @(t) (t>t_glut);
theta2 = @(V,t) (V>-30).*(t>t_bpap).*(t<t_bpap+2);
H=@(V) (V_Ca-V)/(1+(MG2/3.57)*exp(-62*V));
I_NMDA = @(V,t) P0*G_NMDA.*theta1(t).*(I_glut_f*exp((t_glut-t)/tau_glut_f)+I_glut_s*exp((t_glut-t)/tau_glut_s)).*H(V);
I_VDCC= @(V,t) G_VDCC.*theta2(V,t).*(V_VDCC-V);
V_bpap = @(t) Vmax*((I_bpap_f*exp(-(t-t_bpap)/tau_bpap_f))+(I_bpap_s*exp(-(t-t_bpap)/tau_bpap_s)));
g=@(t,Ca) I_NMDA(V,t)+I_VDCC(V,t)-(Ca-Ca_basal)/tau_decay;
[t,Ca] =ode45(g,span,0.1);
plot(t, Ca, 'LineWidth', 1.5);
xlim([1000 2000]);