source "./geometry_helpers_3d.m"
source "./bawm_indices.m"
source "./bawm_landmarks.m"

function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)=pertLand(XL(:,landmark_index),dxl);
  endfor;
endfunction;

function [XR, XL,chi_stats_l, num_inliers_l,H, b]=doTotalLS(XR, XL,
	                                                   Zl,
                                                     landmark_associations,
	                                                   num_poses,
	                                                   num_landmarks,
	                                                   num_iterations,
	                                                   damping,
	                                                   kernel_threshold)

  global pose_dim;
  global landmark_dim;
  global matchable_dim;

  chi_stats_l=zeros(1,num_iterations);
  num_inliers_l=zeros(1,num_iterations);
  chi_stats_p=zeros(1,num_iterations);
  num_inliers_p=zeros(1,num_iterations);
  chi_stats_r=zeros(1,num_iterations);
  num_inliers_r=zeros(1,num_iterations);

  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
   
    if (num_landmarks) 
      [H_landmarks, b_landmarks, chi_, num_inliers_] = linearizeLandmarks(XR, XL, Zl, landmark_associations,num_poses, num_landmarks, kernel_threshold);
      chi_stats_l(iteration)=chi_;
      num_inliers_l(iteration)=num_inliers_;
    endif;
    
    ## H=H_poses;
    ## b=b_poses;
    
    if (num_landmarks) 
       H+=H_landmarks;
       b+=b_landmarks;
    endif;

    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    r = rank(H(pose_dim+1:end,pose_dim+1:end))
    s = size(H(pose_dim+1:end,pose_dim+1:end))
    
    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor

endfunction

## function P_world=makeWorld(num_landmarks,world_size)
##   global matchable_dim;
##   P_world=zeros(matchable_dim,num_landmarks);

##   for i=1:num_landmarks
##     L=zeros(matchable_dim,1);
##     L(1:3,1)     = (rand(3,1)-0.5)*world_size;
##     #Rl = eye(3);
##     dl=rand(3,1);
##     dl/=norm(dl);
##     Rl=d2R(dl);
##     L(4:12,1)    = Rl(:);
##     L(13) = 2;
##     P_world(:,i) = L;
##   endfor
## endfunction

function P_world=makePlanesWorld(num_landmarks,world_size)
  global matchable_dim;
  P_world=zeros(matchable_dim,num_landmarks);

  for i=1:num_landmarks
    L=zeros(matchable_dim,1);
    L(1:3,1)     = (rand(3,1)-0.5)*world_size;
    #Rl = eye(3);
    dl=rand(3,1);
    dl/=norm(dl);
    Rl=d2R(dl);
    L(4:12,1)    = Rl(:);
    L(13) = 3;
    P_world(:,i) = L;
  endfor
endfunction

function [Zl,landmark_associations]=generateMeasurements(num_poses,num_landmarks,XR_true,XL_true)
  global matchable_dim;
  num_landmark_measurements=num_poses*num_landmarks;
  Zl=zeros(matchable_dim,num_landmark_measurements);
  landmark_associations=zeros(2,num_landmark_measurements);
  
  measurement_num=1;
  for pose_num=1:num_poses
    Xr=inv(XR_true(:,:,pose_num));
    for (landmark_num=1:num_landmarks)
	    Xl=XL_true(:,landmark_num);
	    landmark_associations(:,measurement_num)=[pose_num,landmark_num]';

      constraint_type=3;
      
      switch(constraint_type)
        case 1
          Zl(:,measurement_num)=generatePointMeasurement(Xr,Xl);
        case 2
          Zl(:,measurement_num)=generateLineMeasurement(Xr,Xl);
        case 3
          Zl(:,measurement_num)=generatePlaneMeasurement(Xr,Xl);
        otherwise
          disp('irrumati');
      endswitch
          
      measurement_num++;
    endfor;
  endfor
endfunction

function Z=generatePlaneMeasurement(Xr,Xl)

  #transform landmark in robot frame
  Z=transLand(Xl,Xr);

endfunction

function Z=generateLineMeasurement(Xr,Xl)

  #plane defined by robot pose
  pr=Xr(1:3,4);
  nr=Xr(1:3,3);
  dr=-nr'*pr;   

  #plane defined by landmark
  pl=Xl(1:3,1);
  #nl=Xl(10:12,1); 
  nl=Xl(4:6,1); 
  dl=-nl'*pl;   

  #line centroid
  A=[nr(1), nr(2);
     nl(1), nl(2)];
  b=[-dr;
     -dl];
  pz=A\b;

  #line direction
  nz=cross(nr,nl);
  nz/=norm(nz);

  #generate measurement
  global matchable_dim;
  Z=zeros(matchable_dim,1);
  Z(1:3,1)=[pz;0];
  Rz=d2R(nz);
  Z(4:12,1)=Rz(:);
  Z(13)=2;
endfunction

function Z=generatePointMeasurement(Xr,Xl)

  #plane defined by robot pose
  pr=Xr(1:3,4);
  nr=Xr(1:3,3);
  dr=-nr'*pr;   

  #plane defined by landmark
  pl=Xl(1:3,1);
  #nl=Xl(10:12,1); 
  nl=Xl(4:6,1); 
  dl=-nl'*pl;   

  #point
  A=[nr(1), nr(2);
     nl(1), nl(2)];
  b=[-dr;
     -dl];
  pz=A\b;

  #generate measurement
  global matchable_dim;
  Z=zeros(matchable_dim,1);
  Z(1:3,1)=[pz;0];
  Rz=eye(3);
  Z(4:12,1)=Rz(:);
  Z(13)=1;
endfunction
