% function [pts_rot, matrix_rot, cntr_out] = rot3d_plan(pts,plan,alpha,cntr);
%
% 3D rotation of points in a plane.
% Rotation is defined by a 3D plane of rotation, the 3D centre of rotation
% in the plane and the angle of rotation.
% The normal to the plane is defined by the cross product of the vectors
% [plan(:,2) - plan(:,1)] and [plan(:,3) - plan(:,1)]. 
%
% Tnputs
%      pts(nb_pts,iXYZ): 3D coordinates of the points
%      plan(iXYZ,no_pt): Three points of the plane
%      alpha           : Rotation angle (radian) (>0 = anticlockwise) in
%                        the plane normal to the direction vector
%      cntr(1,iXYZ)    : Rotation centre (in the plane)
%
% outputs
%      pts_rot(nb_pts,iXYZ): 3D coordinates of the rotated points
%      matrix_rot(3,3)     : Rotation matrix
%      cntr_out(iXYZ)      : Rotation centre (in the plane) (= cntr)
%
% Author: Antoine Serrurier
% Date (adaptation): 24/06/2024

function [pts_rot, matrix_rot, cntr_out] = rot3d_plan(pts,plan,alpha,cntr)

% Maths: a rotation of matrix R and centre C correspo,ds to the following
% transformation:
%  pts_out = R * (pts_in - C) + C
%
% Here, the rotation is defined by a plane (de normal vector u = (ux, uy,
% uz)), a rotation alpha and a rotation centre in the plane. The rotation
% matrix is defined as follows:
%ct = cos(alpha);
%st = sin(alpha);
%R = [ct+(1-ct)*ux^2      (1-ct)*uy*ux-st*uz  (1-ct)*ux*uz+st*uy;...
%     (1-ct)*ux*uy+st*uz  ct+(1-ct)*uy^2      (1-ct)*uz*uy-st*ux;...
%     (1-ct)*ux*uz-st*uy  (1-ct)*uy*uz+st*ux  ct+(1-ct)*uz^2];


% Centring
pts = pts - repmat(cntr,size(pts,1),1);

% Normal vector
vect1 = plan(:,2) - plan(:,1);
vect2 = plan(:,3) - plan(:,1);
direction = cross(vect1, vect2);
direction = direction./norm(direction);

% Rotation matrix
ct = cos(alpha);
st = sin(alpha);
ux = direction(1);
uy = direction(2);
uz = direction(3);
M = [ct+(1-ct)*ux^2      (1-ct)*uy*ux-st*uz  (1-ct)*ux*uz+st*uy;...
     (1-ct)*ux*uy+st*uz  ct+(1-ct)*uy^2      (1-ct)*uz*uy-st*ux;...
     (1-ct)*ux*uz-st*uy  (1-ct)*uy*uz+st*ux  ct+(1-ct)*uz^2];

% Rotation
pts_rot = M * pts';
pts_rot = pts_rot';

% De-centring
pts_rot = pts_rot + repmat(cntr,size(pts,1),1);


% Output
matrix_rot = M;
cntr_out = cntr;


