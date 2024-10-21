clc
clear
close 
%Initial and Boundary condition
v1=102;
v2=118;
v3=167;
v4=235;
v5=289;
v6=302;
v7=299;
v8=288;
v9=250;
v10=234;
v11=203;
v12=167;
v13=144;
v14=117;
v15=104;
v16=102;
v17=102;
v18=102;
v19=102;
v20=102;
v21=102;
L=30000;
g=9.81; 
n=0.023;
delx=1000;
delt=100;
gama=delt/delx;
kmax=721;
tmax=72000;
N=(L/delx)+1;
i=1:N;
x=0:delx:L;
t=0:delt:tmax;
sb=0.00061;
for k=1:37
    D(1,k)=v1+((v2-v1)/36)*(k-1); 
end
for k=37:73
    D(1,k)=v2+((v3-v2)/36)*(k-37);
end
for k=73:109
    D(1,k)=v3+((v4-v3)/36)*(k-73);
end
for k=109:145
    D(1,k)=v4+((v5-v4)/36)*(k-109);
end
for k=145:181
    D(1,k)=v5+((v6-v5)/36)*(k-145);
end
for k=181:217
    D(1,k)=v6+((v7-v6)/36)*(k-181);
end
for k=217:253
    D(1,k)=v7+((v8-v7)/36)*(k-217);
end
for k=253:289
    D(1,k)=v8 + ((v9-v8)/36)*(k-253);
end
for k = 289:325
    D(1,k) = v9 + ((v10-v9)/36)*(k-289);
end
for k = 325:361
    D(1,k) = v10 + ((v11-v10)/36)*(k-325);
end
for k = 361:397
    D(1,k) = v11 + ((v12-v11)/36)*(k-361);
end
for k = 397:433
    D(1,k) = v12 + ((v13-v12)/36)*(k-397);
end
for k = 433:469
    D(1,k) = v13 + ((v14-v13)/36)*(k-433);
end
for k = 469:505
    D(1,k) = v14 + ((v15-v14)/36)*(k-469);
end
for k = 505:541
    D(1,k) = v15 + ((v16-v15)/36)*(k-505);
end
for k = 541:577
    D(1,k) = v16 + ((v17-v16)/36)*(k-541);
end
for k = 577:613
    D(1,k) = v17 + ((v18-v17)/36)*(k-577);
end
for k = 613:649
    D(1,k) = v18 + ((v19-v18)/36)*(k-613);
end
for k = 649:685
    D(1,k) = v19 + ((v20-v19)/36)*(k-649);
end
for k=685:721
  D(1,k)=v20;
end
b(1,1)=120;

for k=1:721
   h(1,k)=(D(1,k)*n/(b(1,1)*sb^(1/2)))^(3/5);
   u(1,k)=D(1,k)/(b(1,1)*h(1,k));
end
for i=1:N 
h(i,1)=h(1,1);
u(i,1)=u(1,1);
end
tic 

for K=1:kmax-1
     % Predictor step (Forward differences)
    for i = 1:N
        sf(i,k) = ((n^2) * u(i,k)^2) / (h(i,k)^(4/3));  
        s2(i,k) = g * h(i,k) * (sb - sf(i,k));           
    end
    for i=2:N-1
        % Continuity equation - predictor step
        hp(i,k+1) = h(i,k) - gama * (h(i,k) * u(i,k) - h(i-1,k) * u(i-1,k));

        % Momentum equation - predictor step
        hpup(i,k+1) = h(i,k) * u(i,k) - gama * (h(i,k) * u(i,k)^2 + 0.5 * g * h(i,k)^2 - h(i-1,k) * u(i-1,k)^2 - 0.5 * g * h(i-1,k)^2) + s2(i,k) * delt;
        up(i,k+1) = hpup(i,k+1) / hp(i,k+1); % Predicted velocity
    end
    
    % Apply boundary conditions for predictor step
    hp(1,k+1) = h(1,k); 
    up(1,k+1) = u(1,k); 
    % Downstream boundary condition 
    hp(N,k+1) = 2 * hp(N-1,k+1) - hp(N-2,k+1);  % Linear extrapolation for water depth
    up(N,k+1) = 2 * up(N-1,k+1) - up(N-2,k+1);  % Linear extrapolation for velocity
    % Corrector step (Rearward differences)
    for i = 2:N-1
       % Continuity equation - corrector step
       hc(i,k+1) = h(i,k) - gama * (hp(i+1,k+1) * up(i+1,k+1) - hp(i,k+1) * up(i,k+1));
       % Momentum equation - corrector step
       hcuc(i,k+1) = h(i,k) * u(i,k) - gama * (hp(i+1,k+1) * up(i+1,k+1)^2 + 0.5 * g * hp(i+1,k+1)^2 - hp(i,k+1) * up(i,k+1)^2 + 0.5 * g * hp(i,k+1)^2) + s2(i,k) * delt;
       uc(i,k+1) = hcuc(i,k+1) / hc(i,k+1); 
    end

    % Apply boundary conditions for corrector step
    h(1,k+1) = hp(1,k+1); 
    u(1,k+1) = up(1,k+1); 
    % Downstream boundary condition - Apply extrapolation
    h(N,k+1) = 2 * hc(N-1,k+1) - hc(N-2,k+1);  % Linear extrapolation for water depth
    u(N,k+1) = 2 * uc(N-1,k+1) - uc(N-2,k+1);  % Linear extrapolation for velocity
end
toc

for i = 1:N
    for k = 1:kmax
        D(i,k) = (1/n) * b(1,1) * h(i,k)^(5/3) * sb^(1/2); 
    end
end

% Plot Hydrograph 
figure;
plot(t, D(1,:), 'k', 'DisplayName', 'Upstream');
hold on;
plot(t, D(N,:), 'r', 'DisplayName', 'Downstream');
xlabel('Time (s)');
ylabel('Discharge (m^3/s)');
title('Hydrograph (Discharge vs Time)');
legend('Upstream', 'Downstream');
grid on;
hold off; 
correct the above code 
