K0=0.5; k1=2; k2=15;k3=1;k4=120;k11=2;k12=15;k13=1;k14=80;Km1=10;Km2=0.3;Km11=10;Km12=1;
Km=4;Ktot=20;Ptot=20;P0=0.5;Ca=0;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;

%pK and P simulation
%pK=X(1),P=X(2)
span = [0 120];

g=@(t,Y,Ca)[k1*((Ktot-Y(1))/(Km1+(Ktot-Y(1))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca^4)*(Ktot-Y(1)))/(Km^4+Ca^4);
    (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca^3)/(Km^3+Ca^3))*(Ptot-Y(2)))];

%basal state modeling
Ca1 = @(t) 0.1*(t<20) + 0.4 * (t>=20)*(t<22)+0.1*(t>=22)*(t<40)+4*(t>=40)*(t<42)+0.1*(t>=42)*(t<60)+2.2*(t>=60)*(t<61.4)+0.1*(t>=61.4)*(t<80)+2.2*(t>=80)*(t<82)+0.1*(t>=82)*(t<100)+2.93*(t>=100)*(t<102)+0.1*(t>=102)*(t<120);

h1 = @(t,y) g(t,y,Ca1(t));
bas= @(t,y) g(t,y,Ca_basal);


[t1,X] =ode45(h1,span,[0 0]);
[t3,B] =ode45(bas,span,[0 0]);


%Normalized EPSP(t) = A(t)/A_{basal}
%AMPAR basal
pK=@(t) interp1(t3, B(:,1), t);
P=@(t) interp1(t3, B(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
syms t3 A_basal
h1=@(t3,A_basal) k21(t3)*(Atot-A_basal)-k22(t3)*A_basal;
[t3,A_basal] =ode45(h1,span,[0.5 0.5]);


% Create a new figure

% A(t)
pK=@(t) interp1(t1, X(:,1), t);
P=@(t) interp1(t1, X(:,2), t);
k21=@(t) c_1*pK(t)+c_3;
k22=@(t) c_2*P(t)+c_4;
syms t1 AX
h2=@(t1,AX) k21(t1)*(Atot-AX)-k22(t1)*AX;
[t1,AX] =ode45(h2,span,[0.5 0.5]);


% Add legend and labels
plot(t1,AX/A_basal, 'LineWidth', 1.5);
ylim([0 2]);
xlabel('Time (t)');
ylabel('A');