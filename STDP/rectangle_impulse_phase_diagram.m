%This is the code making phase diagram for rectangle impulse: constant calcium value with certain duration. Result is saved in the txt file

a=0.25;

b=0.9 ;

c=0.5;

ratio=100;
Ktot=20;Ptot=20;P0=0.5;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;
K0=c; k1=2/10; k2=15/10; k3=1/10; k4=a*120/10; k11=0.075; k12=15/10; k13=1/10; k14=a*90/10;
Km1=10;Km2=0.3;Km11=Ptot/ratio;Km12=b;
Km=2;
Km4=4;

v=0.5; span = [0 200];
initial_condition=[0.018,0.085,5.708];
fid = fopen('Km=2_k14=2.25.txt', 'w');
%mu_par=[init,threshold,height], we use step function for the mus
%[index,thr,h]
mu_par = [1.2,2,40];
%par2=[beta,lambda]
par2=[0.9,1/5];
for dt=0.1:0.1:10
    %fprintf(fid, 'dt=%.2f, ', dt);
    for impulse=0.1:0.1:10
        result_array = zeros(1, 100);
        Ca = @(t) 0.1 + impulse * (t>8 & t<8+dt);
        mu = @(t) (Ca(t)^mu_par(1))./((Ca(t)^mu_par(1)) + mu_par(2)^mu_par(1))*mu_par(3);
        nu = @(pK) v*par2(1)/(1+ par2(2)*pK)+1-par2(1);
        g = @(t,Y) [k1*((Ktot-Y(1)-Y(3))/(Km1+(Ktot-Y(1)-Y(3))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca(t).^4)*(Ktot-Y(1)-Y(3)))/(Km4^4+Ca(t).^4);
            (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca(t).^3)/(Km^3+Ca(t).^3))*(Ptot-Y(2)));
            -mu(t)*Y(3)+nu(Y(1))*(Ktot-Y(1)-Y(3))];
               
        % Solve ODE
        opts = odeset('RelTol',1e-6,'AbsTol',1e-5, 'MaxStep', 0.1);
        [t, X] = ode45(g, span, initial_condition, opts);


        difference = X(end,1) - X(end,2);
        idx=round(impulse/0.1);
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