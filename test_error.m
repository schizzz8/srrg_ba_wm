close all
clear
clc

source "./bawm_landmarks.m"

#robot pose
#Xr = eye(4);
Xr = v2t([1;2;3;0.1;0.2;0.3]);

## Z = zeros(15,1);
## Z(1:3,1) = [10;5;4];
## #RZ = eye(3);
## RZ = Rx(0.5)*Ry(0.6)*Rz(0.7);
## Z(4:12,1) = RZ(:);
## Z(13:15,1)= [1;1;1];

## Xl = transLand(Z, Xr)

## [e,Jrn,Jln]=landmarkErrorAndJacobian(Xr,Xl,Z)

Z  = [1.31829;-1.85795;-1.95754;1;0;0;0;1;0;0;0;1;1;1;1];

Xl = [1.57544;-1.73499;-1.77925;1;0;0;0;1;0;0;0;1;1;1;1];

R=Xr(1:3,1:3);
t=Xr(1:3,4);

pl=Xl(1:3,1);
Rl=reshape(Xl(4:12,1),3,3);

pz=Z(1:3,1);
RZ=reshape(Z(4:12,1),3,3);

ep = Rl'*(R*pz + t - pl)
ed = (R*RZ - Rl)(:,1)
eo = (R*RZ)(:,1)' * Rl(:,1)

e = [ep;ed;eo];
