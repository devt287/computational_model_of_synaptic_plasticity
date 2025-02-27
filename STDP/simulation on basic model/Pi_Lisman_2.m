K0=0.5; k1=2; k2=15;k3=1;k4=120;k11=2;k12=15;k13=1;k14=80;Km1=10;Km2=0.3;Km11=10;Km12=1;
Km=4;Ktot=20;Ptot=20;P0=0.5;Ca=0;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;

%pK and P simulation
%pK=X(1),P=X(2)
span = [0 20];

g=@(t,Y,Ca)[k1*((Ktot-Y(1))/(Km1+(Ktot-Y(1))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca^4)*(Ktot-Y(1)))/(Km^4+Ca^4);
    (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca^3)/(Km^3+Ca^3))*(Ptot-Y(2)))];

%basal state modeling
Ca1 = @(t) 0.1 + 4 * (t>8)*(t<10);
Ca2=  @(t) 0.1 + 2.2 * (t>8)*(t<10);
h1 = @(t,y) g(t,y,Ca1(t));
h2 = @(t,y) g(t,y,Ca2(t));
bas= @(t,y) g(t,y,Ca_basal);


[t1,X] =ode45(h1,span,[0 0]);
[t2,Y] =ode45(h2,span,[0 0]);
[t3,B] =ode45(bas,span,[0 0]);


%Normalized EPSP(t) = A(t)/A_{basal}
%AMPAR basal
pK=@(t) interp1(t3, B(:,1), t);
P=@(t) interp1(t3, B(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;

h1=@(t3,A_basal) k21(t3)*(Atot-A_basal)-k22(t3)*A_basal;
[t3,A_basal] =ode45(h1,span,[0.5 0.5]);


% Create a new figure

% A(t) LTP
pK=@(t) interp1(t1, X(:,1), t);
P=@(t) interp1(t1, X(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h2=@(t1,AX) k21(t1)*(Atot-AX)-k22(t1)*AX;
[t1,AX] =ode45(h2,span,[0.5 0.5]);


% A(t) LTD
pK=@(t) interp1(t2, Y(:,1), t);
P=@(t) interp1(t2, Y(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
h3=@(t2,AY) k21(t2)*(Atot-AY)-k22(t2)*AY;
[t2,AY] =ode45(h3,span,[0.5 0.5]);



% Add legend and labels
plot(t1,AX/A_basal,t2,AY/A_basal, 'LineWidth', 1.5);
legend('A(t) LTP Ca=4.0', 'A(t) LTD Ca=2.2');
xlabel('Time (t)');
ylabel('A');
ylim([0 2]);












