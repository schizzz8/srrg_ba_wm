#(minimal) size of pose and landmarks
global pose_dim=6;
global landmark_dim=5;
global matchable_dim=13; #last field is the type: 1:point, 2:line, 3:plane

function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  if (pose_index>num_poses)
    v_idx=-1;
    return;
  end
  v_idx=1+(pose_index-1)*pose_dim;
end

function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  end
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
end
