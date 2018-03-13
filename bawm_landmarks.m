source "./geometry_helpers_3d.m"
source "./bawm_indices.m"

function X = landmarkPerturbation(Xl,dx)
  dp = dx(1:3);
  dR = Ry(dx(4))*Rz(dx(5));

  X = eye(4);
  X(1:3,4)   = Xl(1:3,4) + dp;
  X(1:3,1:3) = Xl(1:3,1:3)*dR;

endfunction

function e = errorFunction(Xr,Xl,Z)
  R=Xr(1:3,1:3);
  t=Xr(1:3,4);
  
  pl=Xl(1:3,4)
  Rl=Xl(1:3,1:3)
  
  pz=Z(1:3,4)
  Rz=Z(1:3,1:3)
  
  ep = Rl'*(R*pz + t - pl);
  ed = (R*Rz - Rl)(:,1);
  eo = (R*Rz)(:,1)' * Rl(:,1);

  e = [ep;ed;eo];

endfunction

function [e,Jr,Jl]=landmarkErrorAndJacobian(Xr,Xl,Z)

  e = errorFunction(Xr,Xl,Z);
  
  epsilon = 1e-3;
  
  Jr = zeros(7,6);
  dxr=zeros(6,1);
  for i=1:6
    dxr(i)=epsilon;
    Xr_plus=v2t(dxr)*Xr;
    
    dxr(i)=-epsilon;
    Xr_minus=v2t(dxr)*Xr;

    Jr(:,i)=errorFunction(Xr_plus,Xl,Z) - errorFunction(Xr_minus,Xl,Z);
    dxr(i)=0;

  endfor

  Jr /= (2*epsilon);

  Jl = zeros(7,5);
  dxl=zeros(5,1);
  for i=1:5
    dxl(i)=epsilon;
    Xl_plus=landmarkPerturbation(Xl,dxl);
    
    dxl(i)=-epsilon;
    Xl_minus=landmarkPerturbation(Xl,dxl);
    
    
    Jr(:,i)=errorFunction(Xr,Xl_plus,Z) - errorFunction(Xr,Xl_minus,Z);
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

  
end
