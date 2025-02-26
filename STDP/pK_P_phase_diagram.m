%This is the code making phase diagram for pK versus P. Result is saved in the txt file




a = 0.25;
b = 0.9;
c=0.5;
pre=1.5;
%[1.45,1.35,1.475,1.5]

K0=c; k1=2/10; k2=15/10; k3=1/10; k4=a*120/10; k11=2/10; k12=15/10; k13=1/10; k14=a*80/10;
Km1=10;Km2=0.3;Km11=10;Km12=b;Km=4;Ktot=20;Ptot=20;P0=0.5;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;
v=0.5; span = [0 200];Ca=0.1; 

fid = fopen('pre=1.5.txt', 'w');
%mu_par=[init,threshold,height], we use step function for the mus
mu_par = [1, 2, 40];
%par2=[beta,lambda]
par2=[0.9,1/5];

pK_init = 0.018;
P_init = 0.085;
C_init = 5.708;
C_tot = (k11*(Ptot - P_init)*P_init) / (Km11 + Ptot - P_init);  


for P_init=0:0.2:20
    %fprintf(fid, 'dt=%.2f, ', dt);
    for pK_init=0:0.2:20
        ratio=100;

        Km11=Ptot/ratio;
        %k11_new=C_tot*(Km11+Ptot-0.085)/((Ptot-0.085)*0.085);
        k11_new=0.075;
        initial_condition=[pK_init, P_init];
        result_array = zeros(1, 100);       
        mu = @(t) (Ca^mu_par(1))./((Ca^mu_par(1)) + mu_par(2)^mu_par(1))*mu_par(3);
        nu = @(pK) v*par2(1)/(1+ par2(2)*pK)+1-par2(1);
        C_expr = nu(pK_init)*(Ktot-pK_init-K0)/(mu(Ca_basal));
        g = @(t,Y) [k1*((Ktot-Y(1)-C_expr)/(Km1+(Ktot-Y(1)-C_expr)))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca.^4)*(Ktot-Y(1)-C_expr))/(Km^4+Ca.^4);
            (k11_new*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(pre*Y(1)+K0)+k13*P0+(k14*(Ca.^3)/(Km^3+Ca.^3))*(Ptot-Y(2)))];

               
        % Solve ODE
        opts = odeset('RelTol',1e-6,'AbsTol',1e-5, 'MaxStep', 0.1);
        [t, X] = ode45(g, span, initial_condition, opts);


        difference = X(end,1) - X(end,2);
        idx=round(pK_init/0.2)+1;
        if difference > 1
            result_array(idx) = 1;
        elseif difference < -1
            result_array(idx) = -1;
        else
            result_array(idx) = 0;
        end
        fprintf(fid,'%d ', result_array(idx));

    end
    fprintf(fid, '\n');
    
end