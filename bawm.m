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

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor

endfunction

function P_world=makeWorld(num_landmarks,world_size)
  global matchable_dim;
  P_world=zeros(matchable_dim,num_landmarks);

  for i=1:num_landmarks
    L=zeros(matchable_dim,1);
    L(1:3,1)     = (rand(3,1)-0.5)*world_size;
    Rl = eye(3);
    L(4:12,1)    = Rl(:);
    L(13:15,1)   = ones(3,1);
    P_world(:,i) = L;
  endfor
endfunction
