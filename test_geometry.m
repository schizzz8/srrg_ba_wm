clear
clc
close all

source "./bawm.m"

#robot pose
Xr=eye(4);
Xr(1:3,1:3)=eye(3);
Xr(1:3,4)  =[1;0;0];

#plane landmark
R=Ry(-1.57);
t=[1;0;0];

Xl=zeros(13);
Xl(1:3) =t;
Xl(4:12)=R(:);
Xl(13)  =3;


#generate measurement
Z=measureLineFromPlane(Xr,Xl)
