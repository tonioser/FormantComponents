function DroiteStruct = droiteMatrix2Struct(DroiteMatrix)

% function DroiteStruct = droiteMatrix2Struct(DroiteMatrix);
% 
% Transform a line ('Droite' in French) from a matrix form to a structure
% form.
% 
% Inputs
%      DroiteMatrix(2, iXYZ) : Line in matrix form:
%                              DroiteMatrix(1, iXYZ) : Point of the line
%                              DroiteMatrix(2, iXYZ) : Director vector
% 
% Outputs
%      DroiteStruct(2, iXYZ) : Line in structure form:
%                              DroiteStruct.pts   : Point of the line
%                              DroiteStruct.V_dir : Director vector
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 22/07/2015

DroiteStruct.pts = DroiteMatrix(1,:);
DroiteStruct.V_dir = DroiteMatrix(2,:);

return



