clear all; clc; close all;

Iterations = 10000;

Ts = 10^-2;
t = 0:Iterations-1;
t = t*Ts;

radius = 15;
l = 45;
teta = pi/2;
vel_max = 2*pi;
phi = 0:(vel_max/Iterations):vel_max;

Rot = rotationMatrix(teta);
A = [.5 .5 0; 0 0 1; 1/(2*l) -1/(2*l) 0];


PosI = zeros(3, Iterations); % X, Y, teta
PosR = zeros(3, Iterations);
VelI = zeros(3, Iterations); % Vx, Vy, omega
VelR = zeros(3, Iterations);

ZGps = zeros(3, Iterations);
ZIns = zeros(4, Iterations);

IMU = zeros(3, Iterations);
IMUI = zeros(3, Iterations);
TetaIMU = zeros(1, Iterations);  % NÃ£o Ã© verdade que o teta da IMU Ã© teta mas Ã© o inicial do sistema (Acredito que faz sentido)
TetaIMU(1,1) = teta;
VelIMU = zeros(2, Iterations);
PosIMU = zeros(2, Iterations);
time_passed = 0;
last_gps = 0;

gpsDeviation = 0.001;
imuDeviation = 0.001;

Q = eye(3)*gpsDeviation*10;
R = eye(3)*imuDeviation*.2;

Pk = zeros(5,5); %Covariance
Xe = zeros(5, Iterations);
Ptrace = zeros(1,Iterations);

Z = zeros(2,Iterations);

for i=2:Iterations
  %%Modelo do robo
  U = [radius*phi(i-1) radius*phi(i-1) 0]';
  VelR(:, i) = A*U;
  VelI(:, i) = inv(Rot)*VelR(:,i);

  AcelR = (VelR(:,i) - VelR(:,i-1))/Ts;
  AcelI = (VelI(:,i) - VelI(:,i-1))/Ts;

  PosR(:,i) = PosR(:, i-1) + VelR(:, i-1)*Ts;
  PosI(:,i) = PosI(:, i-1) + VelI(:, i-1)*Ts;

  time_passed = time_passed + Ts;
  if((time_passed - last_gps)>(Ts*100))
      ZGps(:,i) = [PosI(1, i) PosI(2, i) sqrt(VelI(1, i)^2 + VelI(2, i)^2)]' + normrnd(0, gpsDeviation);
  else
      ZGps(:,i) = ZGps(:,i-1);
  end
  

  %INS
  IMU(:,i) = [AcelR(1) AcelR(2) VelR(3, i)]' + normrnd(0, imuDeviation);
  TetaIMU(1, i) = TetaIMU(1, i-1) + IMU(3,i-1)*Ts;
  RIMU = rotationMatrix(TetaIMU(1, i));
  IMUI(:,i) = inv(RIMU)*IMU(:,i);

  VelIMU(:, i) = VelIMU(:, i-1) + IMUI(1:2, i-1)*Ts;
  PosIMU(:, i) = PosIMU(:, i-1) + VelIMU(:, i-1)*Ts;

  ZIns(:,i) = [PosIMU(1, i) PosIMU(2, i) VelIMU(1, i) VelIMU(2, i)]';

  %EFK
  if((time_passed - last_gps)>(Ts*100))
    last_gps = time_passed;
    F = [0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0; 0 0 -IMUI(2, i) 0 0; 0 0 IMUI(1, i) 0 0]; %Conferir sinais e ver se usa i ou i-1
    G = [0 0 0; 0 0 0; 1 0 0; 0 cos(TetaIMU(1,i)) -sin(TetaIMU(1,i)); 0 sin(TetaIMU(1,i)) cos(TetaIMU(1,i))];

    v_sqrt = sqrt(ZIns(3,i)^2 + ZIns(4,1)^2); %ZIns?
    if(v_sqrt==0)
        H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0];
    else
        H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 ZIns(3,i)/v_sqrt ZIns(4,i)/v_sqrt];
    end
    
    D = zeros(3,3);

    [Fd, Gd, Hd, Vd] = c2dm(F,G,H,D,Ts*100);
    
    dZ = [ZGps(1,i)-ZIns(1,i), ZGps(2,i)-ZIns(2,i), ZGps(3,i)-sqrt(ZIns(3,i)^2 +ZIns(4,i)^2)]';
    
    %Pred
    Xk =0;
    Pk = Fd*Pk*Fd' + Gd*Q*Gd';
    
    %Filter update
    K = Pk * H' *inv(H*Pk*H' + R');
    Xe(:,i) = K*dZ;
    Pk = (eye(5) - K * H) * Pk;
    Ptrace(i) = trace(Pk);
  else
    Xe(:,i) = Xe(:, i-1); % 0?
    Ptrace(i) = trace(Pk);

  end
  Z(1,i) = ZIns(1,i) + Xe(1,i);
  Z(2,i) = ZIns(2,i) + Xe(2,i);
end

figure(1)
subplot(2,2,1);
plot(t, ZIns(1,:), ':k', t, ZGps(1,:), 'b',t, Z(1,:), 'r', t, PosI(1,:),'--b');
legend('Ins','GPS','EKF','Real');
xlabel('Time [s]');
ylabel('X');
title('X vs Time');

% figure(2)
subplot(2,2,2);
plot(t, ZIns(2,:), ':k', t, ZGps(2,:), 'b',t, Z(2,:), 'r', t, PosI(2,:),'--b');
legend('Ins','GPS','EKF','Real');
xlabel('Time [s]');
ylabel('Y');
title('Y vs Time');

% figure(3)
subplot(2,2,3);
plot(ZIns(1,:), ZIns(2,:),':k', ZGps(1,:), ZGps(2,:),'b', Z(1,:), Z(2,:),'r', PosI(1,:), PosI(2,:),'--b');
legend('Ins','GPS','EKF','Real')';
xlabel('X');
ylabel('Y');
title('Y vs X');

% figure(4)
subplot(2,2,4);
plot(t,Ptrace);
xlabel('t');
ylabel('Trace P');
title('P vs Time');

for i=1:20:Iterations
    figure(2);
    %plot da trajetória feita pelo robô
    plot(PosI(1,1:i),PosI(2,1:i),'b--','lineWidth',1)

    hold on
    %plot da posição atual do veículo
    plot(PosI(1,i),PosI(2,i),'bo','lineWidth',3)
    %plot da posição atual do veículo
    plot(ZGps(1, 1:i),ZGps(2, 1:i),'ko','MarkerSize',10)
    %plot da trajetória do EKF
    plot(Z(1,1:i),Z(2,1:i),'--r','lineWidth',3)
    %plot da trajetória do dead reckoning
    plot(ZIns(1,1:i),ZIns(2,1:i),':k','lineWidth',3)
    
    hold off
    xlim([PosI(1,i)-5 PosI(1,i)+5])
    ylim([PosI(2,i)-5 PosI(2,i)+5])
    xlabel('x(m)','FontSize',14) 
    ylabel('y(m)','FontSize',14)
    set(gca, 'FontSize',18)
    drawnow
    pause(.0000000001)

    
end