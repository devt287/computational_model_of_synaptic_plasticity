K0=0.5; k1=2; k2=15;k3=1;k4=120;k11=2;k12=15;k13=1;k14=80;Km1=10;Km2=0.3;Km11=10;Km12=1;
Km=4;Ktot=20;Ptot=20;P0=0.5;Ca=0;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;

%pK and P simulation
%pK=X(1),P=X(2)
g=@(t,Y,Ca)[k1*((Ktot-Y(1))/(Km1+(Ktot-Y(1))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca^4)*(Ktot-Y(1)))/(Km^4+Ca^4);
    (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca^3)/(Km^3+Ca^3))*(Ptot-Y(2)))];

Ca1 = @(t) 0.1 + 0.1 * (t>8)*(t<10);
Ca2=  @(t) 0.1 + 0.2 * (t>8)*(t<10);
Ca3=  @(t) 0.1 + 0.3 * (t>8)*(t<10);
Ca4=  @(t) 0.1 + 0.4 * (t>8)*(t<10);
Ca5=  @(t) 0.1 + 0.5 * (t>8)*(t<10);
Ca6=  @(t) 0.1 + 0.6 * (t>8)*(t<10);

h1 = @(t,y) g(t,y,Ca1(t));
h2 = @(t,y) g(t,y,Ca2(t));
h3 = @(t,y) g(t,y,Ca3(t));
h4 = @(t,y) g(t,y,Ca4(t));
h5 = @(t,y) g(t,y,Ca5(t));
h6 = @(t,y) g(t,y,Ca6(t));
bas= @(t,y) g(t,y,Ca_basal);

span = [0 20]; 

[t0,B] =ode45(bas,span,[0 0]);
[t1,X] =ode45(h1,span,[0 0]);
[t2,Y] =ode45(h2,span,[0 0]);
[t3,Z] =ode45(h3,span,[0 0]);
[t4,N] =ode45(h4,span,[0 0]);
[t5,G] =ode45(h5,span,[0 0]);
[t6,Q] =ode45(h6,span,[0 0]);



%Normalized EPSP(t) = A(t)/A_{basal}

%AMPAR basal
pK=@(t) interp1(t0, B(:,1), t);
P=@(t) interp1(t0, B(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h1=@(t0,A_basal) k21(t0)*(Atot-A_basal)-k22(t0)*A_basal;
[t0,A_basal] =ode45(h1,span,[0.5 0.5]);


% Create a new figure

% perturbation no greater than Ca=0.5, stable
pK=@(t) interp1(t1, X(:,1), t);
P=@(t) interp1(t1, X(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h1=@(t1,AX) k21(t1)*(Atot-AX)-k22(t1)*AX;
[t1,AX] =ode45(h1,span,[0.5 0.5]);

pK=@(t) interp1(t2, Y(:,1), t);
P=@(t) interp1(t2, Y(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h2=@(t2,AY) k21(t2)*(Atot-AY)-k22(t2)*AY;
[t2,AY] =ode45(h2,span,[0.5 0.5]);

pK=@(t) interp1(t3, Z(:,1), t);
P=@(t) interp1(t3, Z(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h3=@(t3,AZ) k21(t3)*(Atot-AZ)-k22(t3)*AZ;
[t3,AZ] =ode45(h3,span,[0.5 0.5]);

pK=@(t) interp1(t4, N(:,1), t);
P=@(t) interp1(t4, N(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h4=@(t4,AN) k21(t4)*(Atot-AN)-k22(t4)*AN;
[t4,AN] =ode45(h4,span,[0.5 0.5]);

pK=@(t) interp1(t5, G(:,1), t);
P=@(t) interp1(t5, G(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h5=@(t5,AG) k21(t5)*(Atot-AG)-k22(t5)*AG;
[t5,AG] =ode45(h5,span,[0.5 0.5]);

% A(t) LTD
pK=@(t) interp1(t6, Q(:,1), t);
P=@(t) interp1(t6, Q(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h6=@(t6,AQ) k21(t6)*(Atot-AQ)-k22(t6)*AQ;
[t6,AQ] =ode45(h6,span,[0.5 0.5]);




% Add legend and labels
hold on;
%plot(t1,AX/A_basal,'LineWidth', 1.5);
plot(t2,AY/A_basal,'LineWidth', 1.5);
plot(t3,AZ/A_basal,'LineWidth', 1.5);
plot(t4,AN/A_basal,'LineWidth', 1.5);
plot(t5,AG/A_basal,'LineWidth', 1.5);
plot(t6,AQ/A_basal,'LineWidth', 1.5);

hold off;
xlabel('Time (t)');
ylabel('A');












