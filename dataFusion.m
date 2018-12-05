clear all; clc; close all;

Iterations = 10000;

Ts = 10^-2;
t = 0:Iterations-1;
t = t*Ts;

gpsDeviation = 0.5;
imuDeviation = 0.9;

radius = 15;
l = 45;
teta = pi/2;
vel_max = 20*pi;
phi = 0:(vel_max/Iterations):vel_max;

R = rotationMatrix(teta);
A = [.5 .5 0; 0 0 1; 1/(2*l) -1/(2*l) 0];


PosI = zeros(3, Iterations); % X, Y, teta
PosR = zeros(3, Iterations);
VelI = zeros(3, Iterations); % Vx, Vy, omega
VelR = zeros(3, Iterations);

ZGps = zeros(3, Iterations);
ZIns = zeros(4, Iterations);

IMU = zeros(3, Iterations);
IMUI = zeros(3, Iterations);
TetaIMU = zeros(1, Iterations);  % Não é verdade que o teta da IMU é teta mas é o inicial do sistema (Acredito que faz sentido)
TetaIMU(1,1) = teta;
VelIMU = zeros(2, Iterations);
PosIMU = zeros(2, Iterations);
time_passed = 0;
last_gps = 0;

P = zeros(5,5); %Covariance
Xe = zeros(5, Iterations);
Ptrace = zeros(1,Iterations);

for i=2:Iterations
  %%Modelo do robo
  U = [radius*phi(i-1) radius*phi(i-1) 0]';
  VelR(:, i) = A*U;
  VelI(:, i) = inv(R)*VelR(:,i);

  AcelR = (VelR(:,i) - VelR(:,i-1))/Ts;
  AcelI = (VelI(:,i) - VelI(:,i-1))/Ts;

  PosR(:,i) = PosR(:, i-1) + VelR(:, i-1)*Ts;
  PosI(:,i) = PosI(:, i-1) + VelI(:, i-1)*Ts;

  time_passed = time_passed + Ts;
  if((time_passed - last_gps)>(Ts*100))
      ZGps(:,i) = [PosI(1, i-1) PosI(2, i-1) sqrt(VelI(1, i-1)^2 + VelI(2, i-1)^2)]' + normrnd(0, gpsDeviation);
      % last_gps = time_passed;
  else
      ZGps(:,i) = ZGps(:,i-1);
  endif

  %INS
  IMU(:,i-1) = [AcelR(1) AcelR(2) VelR(3, i-1)]' + normrnd(0, imuDeviation);
  TetaIMU(1, i) = TetaIMU(1, i-1) + IMU(3,i-1)*Ts;
  RIMU = rotationMatrix(TetaIMU(1, i-1));
  IMUI(:,i-1) = inv(RIMU)*IMU(:,i-1);

  VelIMU(:, i) = VelIMU(:, i-1) + IMUI(1:2, i-1)*Ts;
  PosIMU(:, i) = PosIMU(:, i-1) + VelIMU(:, i-1)*Ts;

  ZIns(:,i) = [PosIMU(1, i-1) PosIMU(2, i-1) VelIMU(1, i-1) VelIMU(2, i-1)]';

  %EFK
  if((time_passed - last_gps)>(Ts*100))
    Xe(:,i) = Xe(:, i-1);
    Ptrace(i) = trace(P);
  else
    F = [0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0; 0 0 -IMUI(2, i-1) 0 0; 0 0 IMUI(1, i-1) 0 0]; %Conferir sinais e ver se usa i ou i-1
    G = [0 0 0; 0 0 0; 1 0 0; 0 cos(TetaIMU(1,i-1)) -sin(TetaIMU(1,i-1)); 0 sin(TetaIMU(1,i-1)) cos(TetaIMU(1,i-1))];

    v_sqrt = sqrt(ZIns(3,i-1)^2 + ZIns(4,i-1)^2); %ZIns?
    if(v_sqrt==0)
        H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0];
    else
        H = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 ZIns(3,i-1)/v_sqrt ZIns(4,i-1)/v_sqrt];
    endif
    D = zeros(3,3);

    [Fd, Gd, Hd, Dd] = c2d(F,G,H,D,Ts*100);
  endif
endfor

plot(t, ZIns(2,:), 'k');
hold
plot(t, ZGps(2,:), 'b');
