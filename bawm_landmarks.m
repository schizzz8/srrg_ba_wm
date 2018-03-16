source "./geometry_helpers_3d.m"
source "./bawm_indices.m"

function Xlnew = pertLand(Xl,dx)
  dp = dx(1:3,1);
  dR = Ry(dx(4))*Rz(dx(5));

  global matchable_dim;
  Xlnew = zeros(matchable_dim,1);
  Xlnew(1:3,1)   = Xl(1:3,1) + dp;
  Xlnew(4:12,1)  = (reshape(Xl(4:12,1),3,3)*dR)(:);
  Xlnew(13:15,1) = Xl(13:15,1);
endfunction

function Xlnew = transLand(Xl,X)
  pl = Xl(1:3,1);
  dl = Xl(4:12,1);
  Rl = reshape(dl,3,3);
  
  t = X(1:3,4);
  R = X(1:3,1:3);

  global matchable_dim;
  Xlnew = zeros(matchable_dim,1);
  Xlnew(1:3,1) = R*pl + t;
  Xlnew(4:12,1) = (R*Rl)(:);
  Xlnew(13:15,1) = Xl(13:15,1);
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

function [e,Jr,Jl]=landmarkErrorAndJacobian(Xr,Xl,Z)

  e = computeError(Xr,Xl,Z);
  
  epsilon = 1e-4;
  
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
end

function [e,Jra,Jla]=landmarkErrorAndAnalyticJacobian(Xr,Xl,Z)

  e = computeError(Xr,Xl,Z);
  
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

function Omega=computeOmega(Xl,Z)
  Omega=zeros(7);

  if(Xl(13:15) == [1;1;1])
    if(Z(13:15) == [1;1;1]) #point-point
      Omega(1:3,1:3) = eye(3);
    endif
    if(Z(13:15) == [0;1;1]) #point-line
      Omega(1:3,1:3) = diag([0;1;1]); 
    endif
    if(Z(13:15) == [1;0;0]) #point-plane
      Omega(1:3,1:3) = diag([1;0;0]);
    endif
  endif

  if(Xl(13:15) == [0;1;1])
    if(Z(13:15) == [1;1;1]) #line-point
      Omega(1:3,1:3) = diag([0;1;1]);
    endif
    if(Z(13:15) == [0;1;1]) #line-line
      Omega(1:3,1:3) = diag([0;1;1]);
      Omega(4:6,4:6) = eye(3)
    endif
    if(Z(13:15) == [1;0;0]) #line-plane
      Omega(1:3,1:3) = diag([1;0;0]);
      Omega(7,7) = 1;
    endif
  endif

  if(Xl(13:15) == [1;0;0])
    if(Z(13:15) == [1;1;1]) #plane-point
      Omega(1:3,1:3) = diag([1;0;0]);
    endif
    if(Z(13:15) == [0;1;1]) #plane-line
      Omega(1:3,1:3) = diag([1;0;0]);
      Omega(7,7) = 1;
    endif
    if(Z(13:15) == [1;0;0]) #plane-plane
      Omega(1:3,1:3) = diag([1;0;0]);
      Omega(4:6,4:6) = eye(3)
    endif
  endif

endfunction

function [H,b, chi_tot, num_inliers]=linearizeLandmarks(XR, XL, Zl, associations,num_poses, num_landmarks, kernel_threshold)

  global pose_dim;
  global landmark_dim;

  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);

  chi_tot=0;
  num_inliers=0;

  for (measurement_num=1:size(Zl,2))
    pose_index=associations(1,measurement_num);
    landmark_index=associations(2,measurement_num);

    z=Zl(:,measurement_num);
    Xr=XR(:,:,pose_index);
    Xl=XL(:,landmark_index);
    [e,Jr,Jl] = landmarkErrorAndAnalyticJacobian(Xr, Xl, z);
    chi=e'*e;
    if (chi>kernel_threshold)
      e*=sqrt(kernel_threshold/chi);
      chi=kernel_threshold;
    else
      num_inliers++;
    endif;
    chi_tot+=chi;
    Omega=computeOmega(Xl,z);

    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*Omega*Jr;

    H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jr'*Omega*Jl;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*Omega*Jl;

    H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jl'*Omega*Jr;

    b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=Jr'*Omega*e;
    b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Jl'*Omega*e;
  endfor
endfunction
