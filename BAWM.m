close all
clear
clc

source "./bawm.m"

# synthesis of the virtual world
num_landmarks=100;
num_poses=10;
world_size=10;

# landmarks in a matrix, one per column
P_world=(rand(landmark_dim, num_landmarks)-0.5)*world_size;

# poses in an array of 4x4 homogeneous transform matrices
XR_true=zeros(4,4,num_poses);
XL_true=P_world;

# initialize 1st pose
XR_true(:,:,1)=eye(4);

# scaling coefficient for uniform random pose generation
# adjusts the translation to cover world_size
# adjusts the rotation to span the three angles;
rand_scale=eye(6);
rand_scale(1:3,1:3)*=(0.5*world_size);
rand_scale(4:6,4:6)*=pi;
for (pose_num=2:num_poses)
    xr=rand(6,1)-0.5;
    Xr=v2t(rand_scale*xr);
    XR_true(:,:,pose_num)=Xr;
endfor;

######################################## LANDMARK MEASUREMENTS ######################################## 
# generate an ideal number of landmark measurements
# each pose observes each landmark
num_landmark_measurements=num_poses*num_landmarks;
Zl=zeros(landmark_dim,num_landmark_measurements);
landmark_associations=zeros(2,num_landmark_measurements);

measurement_num=1;
for (pose_num=1:num_poses)
    Xr=inv(XR_true(:,:,pose_num));
    for (landmark_num=1:num_landmarks)
	    Xl=XL_true(:,landmark_num);
	    landmark_associations(:,measurement_num)=[pose_num,landmark_num]';
      Zl(:,measurement_num)=transLand(Zl(:,measurement_num),Xr);
      measurement_num++;
    endfor;
endfor

############################## GENERATION OF (WRONG) INITIAL GUESS ################################## 

# apply a perturbation to each ideal pose (construct the estimation problem)
pert_deviation=1;
pert_scale=eye(6)*pert_deviation;
XR_guess=XR_true;
XL_guess=XL_true;

for (pose_num=2:num_poses)
    xr=rand(6,1)-0.5;
    dXr=v2t(pert_scale*xr);
    XR_guess(:,:,pose_num)=dXr*XR_guess(:,:,pose_num);
endfor;

#apply a perturbation to each landmark
for landmark_num=1:num_landmarks
  dXl=(rand(5,1)-0.5)*pert_deviation;
  XL_guess(:,landmark_num)=pertLand(XL_true(:,landmark_num),dXl);
endfor

