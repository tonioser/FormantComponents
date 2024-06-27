function circle = circle_least_square(M,varargin)
% 
% Calculate the circle fitting a set of points M in the least-square sense.
%
% Inputs
%     M(nbPts,2) : Coordinates of the 2D points
%   
% Outputs
%     circle
%       circle.centre(1,2): Coordinates of the centre of the circle
%       circle.radius     : Radius of the circle
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 21/06/2024

% Number of points
N = size(M,1);

% Initial solution: least-square on the implicite function of the sphere
Mat = [2*M,ones(N,1)] ;      % Calculation matrix
Vec = sqrt(sum(M.^2,2)).^2 ; % Target vector
% Calculation
result = inv(Mat'*Mat)*Mat'*Vec ;
circle.centre = result(1:2)' ; % Centre initial
circle.radius = sqrt(result(3)+norm(circle.centre)^2) ; % Rayon initial

% Optimisation by linear least square (on distance points-circle)
result = 1;
while norm(result) > 1e-10
    % Formatting
    OMi = M - ones(N,1) * circle.centre ; % Centred points
    nOMi = sqrt(sum(OMi.^2,2)) ;          % Distance (norm)
    Ui = 1 - circle.radius ./ nOMi ;      % Relative differece radius-nOMi
    % Calculation matrix
    L = sum( OMi ./ (nOMi * [1,1])) ;
    Mat = [N,L;L',N*eye(2)] ;
    % Calculation vector
    Vec = [sum(Ui.*nOMi),sum((Ui * [1,1]) .* OMi)]' ;
    % Result
    result = inv(Mat)*Vec ;
    % New centre and radius
    circle.radius = circle.radius + result(1) ;
    circle.centre = circle.centre + result(2:3)' ;
end
