close all
clear
clc

source "./bawm_landmarks.m"

#robot pose
#Xr = eye(4);
Xr = v2t([1;2;3;0.1;0.2;0.3]);


Z = zeros(15);
Z(1:3,1) = [10;5;4];
#RZ = eye(3);
RZ = Rx(0.5)*Ry(0.6)*Rz(0.7);
Z(4:12,1) = RZ(:);

Xl = transLand(Z, Xr);



[e,Jrn,Jln]=landmarkErrorAndJacobian(Xr,Xl,Z)
