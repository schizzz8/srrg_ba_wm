close all
clear
clc

source "./bawm_landmarks.m"

Xr = eye(4);

Xf = eye(4);
Xf(1:3,4) = [5;5;5];
Xf(1:3,1:3) = diag([1,2,3]);

Xm = eye(4);
Xm(1:3,4) = [5;5;5];
Xm(1:3,1:3) = diag([1,2,3]);

[e,Jr,Jl]=landmarkErrorAndJacobian(Xr,Xf,Xm);
