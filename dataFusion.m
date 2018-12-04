clear all; clc; close all;

Iterations = 10000;

Ts = 10^-2;
t = 0:Iterations-1;
t = t*Ts;

gpsDeviation = 0.5;
imuDeviation = 0.5;

radius = 15;
l = 45;
teta = pi/2;
vel_max = 10*pi;
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

for i=2:Iterations
  %%Modelo do robo
  U = [radius*phi(i-1) radius*phi(i-1) 0]';
  VelR(:, i) = A*U;
  VelI(:, i) = inv(R)*VelR(:,i);

  AcelR = (VelR(:,i) - VelR(:,i-1))/Ts;
  AcelI = (VelI(:,i) - VelI(:,i-1))/Ts;

  PosR(:,i) = PosR(:, i-1) + VelR(:, i-1)*Ts;
  PosI(:,i) = PosI(:, i-1) + VelI(:, i-1)*Ts;

  %GPS
  ZGps(:,i) = [PosI(1, i-1) PosI(2, i-1) sqrt(VelI(1, i-1)^2 + VelI(2, i-1)^2)]' + normrnd(0, gpsDeviation);

  %INS
  IMU(:,i-1) = [AcelR(1) AcelR(2) VelR(3, i-1)]' + normrnd(0, imuDeviation);
  TetaIMU(1, i) = TetaIMU(1, i-1) + IMU(3,i-1)*Ts;
  RIMU = rotationMatrix(TetaIMU(1, i-1));
  IMUI(:,i-1) = inv(RIMU)*IMU(:,i-1);

  VelIMU(:, i) = VelIMU(:, i-1) + IMUI(1:2, i-1)*Ts;
  PosIMU(:, i) = PosIMU(:, i-1) + VelIMU(:, i-1)*Ts;

  ZIns(:,i) = [PosIMU(1, i-1) PosIMU(2, i-1) VelIMU(1, i-1) VelIMU(2, i-1)]';

endfor

plot(t, ZIns(2,:), 'k');
hold
plot(t, ZGps(2,:), 'b');
