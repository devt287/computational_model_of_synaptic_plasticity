K0=0.5; k1=2; k2=15;k3=1;k4=120;k11=2;k12=15;k13=1;k14=80;Km1=10;Km2=0.3;Km11=10;Km12=1;
Km=4;Ktot=20;Ptot=20;P0=0.5;Ca=0;Ca_basal=0.1;Atot=1;c_1=1;c_2=1;c_3=6;c_4=8;



%pK and P simulation
%pK=X(1),P=X(2)
span = [0 20];

g=@(t,Y,Ca)[k1*((Ktot-Y(1))/(Km1+(Ktot-Y(1))))*Y(1)-((k2*Y(1))/(Km2+Y(1)))*(Y(2)+P0)+k3*K0+(k4*(Ca^4)*(Ktot-Y(1)))/(Km^4+Ca^4);
    (k11*((Ptot-Y(2))/(Km11+(Ptot-Y(2))))*Y(2)-k12*(Y(2)/(Km12+Y(2)))*(Y(1)+K0)+k13*P0+(k14*(Ca^3)/(Km^3+Ca^3))*(Ptot-Y(2)))];
Ca1 = @(t) 0.1 + 4 * (t>8)*(t<10);
Ca2=  @(t) 0.1 + 2.2 * (t>8)*(t<10);
h1 = @(t,y) g(t,y,Ca1(t));
h2 = @(t,y) g(t,y,Ca2(t));
[t1,X]=ode45(h1,span,[0 0]);
[t2,Y]=ode45(h2,span,[0 0]);


% Create a new figure
figure

% Plot for LTP
subplot(2,2,1); 
hold on;
plot(t1, X(:,1), 'LineWidth', 1.5,'DisplayName','K');
plot(t1, X(:,2), 'LineWidth', 1.5,'DisplayName','P');
title('LTP');
legend('pK','P');
xlabel('Time (t)');
ylabel('pK');
xlim([0 20]);
ylim([-2 20]);

% Plot for LTD
subplot(2,2,2); 
hold on;
plot(t2, Y(:,1), 'LineWidth', 1.5,'DisplayName','K');
plot(t2, Y(:,2), 'LineWidth', 1.5,'DisplayName','P');
title('LTD');
legend('pK','P');
xlabel('Time (t)');
ylabel('pK');
xlim([0 20]);
ylim([-2 20]);