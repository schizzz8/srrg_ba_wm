source "./geometry_helpers_3d.m"
source "./bawm_indices.m"

function Xlnew = pertLand(Xl,dx)
  dp = dx(1:3);
  dR = Ry(dx(4))*Rz(dx(5));

  Xlnew = zeros(9);
  Xlnew(1:3) = Xl(1:3) + dp;
  Xlnew(4:6) = dR*Xl(4:6);
  Xlnew(7:9) = Xl(7:9);
endfunction

function Xlnew = transLand(Xl,X)
  p = Xl(1:3);
  d = Xl(4:6);
  
  t = X(1:3,4);
  R = X(1:3,1:3);

  Xlnew = zeros(9);
  Xlnew(1:3) = R*p + t;
  Xlnew(4:6) = R*d;
  Xlnew(7:9) = Xl(7:9);
endfunction

function e = computeError(Xr,Xl,Z)
  R=Xr(1:3,1:3);
  t=Xr(1:3,4);
  
  pl=Xl(1:3);
  dl=Xl(4:6);
  Rl=d2R(dl);
  
  pz=Z(1:3);
  dz=Z(4:6);
  Rz=d2R(dz);
  
  ep = Rl'*(R*pz + t - pl);
  ed = (R*Rz - Rl)(:,1);
  eo = (R*Rz)(:,1)' * Rl(:,1);

  e = [ep;ed;eo];

endfunction

function [e,Jr,Jl]=landmarkErrorAndJacobian(Xr,Xl,Z)

  e = computeError(Xr,Xl,Z);
  
  epsilon = 1e-3;
  
  Jr = zeros(7,6);
  dxr=zeros(6,1);
  for i=1:6
    dxr(i)=epsilon;
    Xr_plus=v2t(dxr)*Xr;
    
    dxr(i)=-epsilon;
    Xr_minus=v2t(dxr)*Xr;

    Jr(:,i)=computeError(Xr_plus,Xl,Z) - computeError(Xr_minus,Xl,Z);
    dxr(i)=0;

  endfor
  Jr /= (2*epsilon);

  Jl = zeros(7,5);
  dxl=zeros(5,1);
  for i=1:5
    dxl(i)=epsilon;
    Xl_plus=pertLand(Xl,dxl);
    
    dxl(i)=-epsilon;
    Xl_minus=pertLand(Xl,dxl);
    
    
    Jl(:,i)=computeError(Xr,Xl_plus,Z) - computeError(Xr,Xl_minus,Z);
    dxr(i)=0;

  endfor
  Jl /= (2*epsilon);

  ## dep_dt  = Rf';
  ## dep_dR  = -Rf'*skew(R*pz + t);
  
  ## dep_dpf = -Rf';
  ## dep_dRf = -clippedSkew(Rf'*(R*pz + t - pf));

  ## ded_dt  = zeros(3);
  ## ded_dR  = -skew((R*Rz)(:,1));

  ## ded_dpf = zeros(3);
  ## ded_dRf = Rf*clippedSkew([1;0;0]);

  ## de0_dt  = zeros(3);
  ## deo_dR  = -[1,0,0]*Rm'*R'*skew(R*Rm[1;0;0]);

  ## deo_dpf = zeros(3);
  ## deo_dRf = -[1,0,0]*Rm'*R'*Rf*clippedSkew([1;0;0]);
end
