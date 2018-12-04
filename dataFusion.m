clear all; clc; close all;

Iterations = 1000;
Ts = 10^-2;

gpsDeviation = 0.05;
imuDeviation = 0.5;

radius = 15;
l = 45;
teta = pi/2;
phi = 2*pi;

R = rotationMatrix(teta);
A = [.5 .5 0; 0 0 1; 1/(2*l) -1/(2*l) 0];
U = [radius*phi radius*phi 0]';

PosI_0 = [0 0 0]'; % X, Y, teta
PosR_0 = [0 0 0]';
VelI_0 = [0 0 0]'; % Vx, Vy, omega
VelR_0 = [0 0 0]';



TetaIMU_0 = teta; % Não é verdade que o teta da IMU é teta mas é o inicial do sistema (Acredito que faz sentido)
VelIMU_0 = [0 0]';
PosIMU_0 = [0 0]';

for i=2:Iterations
  %%Modelo do robo
  VelR_1 = A*U;
  VelI_1 = inv(R)*VelR_1;

  AcelR = (VelR_1 - VelR_0)/Ts;
  AcelI = (VelI_1 - VelI_0)/Ts;

  PosR_1 = PosI_0 + VelR_0*Ts;
  PosI_1 = PosI_0 + VelI_0*Ts;

  %GPS
  ZGps = [PosI_0(1) PosI_0(2) sqrt(VelI_0(1)^2 + VelI_0(2)^2)]' + normrnd(0, gpsDeviation);

  %INS
  IMU = [AcelR(1) AcelR(2) VelR_0(3)]' + normrnd(0, imuDeviation);
  TetaIMU_1 = TetaIMU_0 + IMU(3)*Ts;
  RIMU = rotationMatrix(TetaIMU_0);
  IMUI = inv(RIMU)*IMU;

  VelIMU_1 = VelIMU_0 + IMUI(1:2)*Ts;
  PosIMU_1 = PosIMU_0 + VelIMU_0*Ts;

  ZIns = [PosIMU_0(1) PosIMU_0(2) VelIMU_0(1) VelIMU_0(2)]';

  %Update Variables
  PosI_0 = PosI_1;
  PosR_0 = PosR_1;

  VelI_0 = VelI_1;
  VelR_0 = VelR_1;

  TetaIMU_0 = TetaIMU_1;
  VelIMU_0 = VelIMU_1;
  PosIMU_0 = PosIMU_1;
endfor
