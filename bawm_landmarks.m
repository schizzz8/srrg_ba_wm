source "./geometry_helpers_3d.m"
source "./bawm_indices.m"

function Xlnew = pertLand(Xl,dx)
  dp = dx(1:3,1);
  dR = Ry(dx(4))*Rz(dx(5));

  Xlnew = zeros(12);
  Xlnew(1:3,1) = Xl(1:3,1) + dp;
  Xlnew(4:12,1) = (reshape(Xl(4:12,1),3,3)*dR)(:);
endfunction

function Xlnew = transLand(Xl,X)
  pl = Xl(1:3,1);
  dl = Xl(4:12,1);
  Rl = reshape(dl,3,3);
  
  t = X(1:3,4);
  R = X(1:3,1:3);

  Xlnew = zeros(12);
  Xlnew(1:3,1) = R*pl + t;
  Xlnew(4:12,1) = (R*Rl)(:);
endfunction

function e = computeError(Xr,Xl,Z)
  R=Xr(1:3,1:3);
  t=Xr(1:3,4);
  
  pl=Xl(1:3,1);
  Rl=reshape(Xl(4:12,1),3,3);
  
  pz=Z(1:3,1);
  RZ=reshape(Z(4:12,1),3,3);
  
  ep = Rl'*(R*pz + t - pl);
  ed = (R*RZ - Rl)(:,1);
  eo = (R*RZ)(:,1)' * Rl(:,1);

  e = [ep;ed;eo];

endfunction

function [e,Jrn,Jln,Jra,Jla]=landmarkErrorAndJacobian(Xr,Xl,Z)

  e = computeError(Xr,Xl,Z);
  
  epsilon = 1e-4;
  
  Jrn = zeros(7,6);
  dxr=zeros(6,1);
  for i=1:6
    dxr(i)=epsilon;
    Xr_plus=v2t(dxr)*Xr;
    
    dxr(i)=-epsilon;
    Xr_minus=v2t(dxr)*Xr;

    Jrn(:,i)=computeError(Xr_plus,Xl,Z) - computeError(Xr_minus,Xl,Z);
    dxr(i)=0;
  endfor
  Jrn /= (2*epsilon);

  Jln = zeros(7,5);
  dxl=zeros(5,1);
  for i=1:5
    dxl(i)=epsilon;
    Xl_plus=pertLand(Xl,dxl);
    
    dxl(i)=-epsilon;
    Xl_minus=pertLand(Xl,dxl);    
    
    Jln(:,i)=computeError(Xr,Xl_plus,Z) - computeError(Xr,Xl_minus,Z);
    dxr(i)=0;
  endfor
  Jln /= (2*epsilon);

  R=Xr(1:3,1:3);
  t=Xr(1:3,4);
  
  pl=Xl(1:3,1);
  Rl=reshape(Xl(4:12,1),3,3);
  
  pz=Z(1:3,1);
  RZ=reshape(Z(4:12,1),3,3);

  dep_dt  = Rl';
  dep_dR  = -Rl'*skew(R*pz + t);
  

  ded_dt  = zeros(3);
  ded_dR  = -skew((R*RZ)(:,1));

  deo_dt  = zeros(1,3);
  deo_dR  = [1,0,0]*RZ'*R'*skew(Rl*[1;0;0]);

  Jra = zeros(7,6);
  Jra = [dep_dt,dep_dR;
         ded_dt,ded_dR;
         deo_dt,deo_dR];

  dep_dpl = -Rl';
  dep_dRl = clippedSkew(Rl'*(R*pz + t - pl));

  ded_dpl = zeros(3);
  ded_dRl = Rl*clippedSkew([1;0;0]);

  deo_dpl = zeros(1,3);
  deo_dRl = -[1,0,0]*RZ'*R'*Rl*clippedSkew([1;0;0]);


  Jla = zeros(7,5);
  Jla = [dep_dpl,dep_dRl;
         ded_dpl,ded_dRl;
         deo_dpl,deo_dRl];
end
