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

function [XR, XL, chi_stats_l, num_inliers_l, chi_stats_p,
          num_inliers_p,chi_stats_r, num_inliers_r, H, b]=doTotalLS(XR, XL,
	                                                                    Zl,
                                                                      landmark_associations,
	     Zp, projection_associations,
	     Zr, pose_associations,
	     num_poses,
	     num_landmarks,
	     num_iterations,
	     damping,
	     kernel_threshold)

  global pose_dim;
  global landmark_dim;

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

    

  endfor


endfunction

