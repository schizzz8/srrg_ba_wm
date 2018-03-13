close all
clear
clc

source "./bawm_landmarks.m"

Xl = eye(4);
Xl(1:3,4) = [0;0;0];
Xl(1:3,1:3) = diag([1,1,1]);

pl = Xl(1:3,4);
Rl = Xl(1:3,1:3);

dx = [0;0;0;1.57;1.57];


dp = dx(1:3);
dR = Ry(dx(4))*Rz(dx(5));

plnew = pl + dp
Rlnew = Rl*dR
