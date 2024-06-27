function Pt_inter = intersection_cnt_cnt_2D(Cnt1, Cnt2)
% 
% function Pt_inter = intersection_cnt_cnt_2D(Cnt1, Cnt2);
% 
% Returns the intersection points between two 2D curves
% 
% Inputs
%   Cnt1(nbPts1, iXY): Coordinates of the points of the first curve
%   Cnt2(nbPts2, iXY): Coordinates of the points of the second curve
% 
% Outputs
%   Pt_inter(nbPtInter, iXY): Coordinates of the intersection points
% 
% Author: Antoine Serrurier
% Date (adaptation): 24/06/2024

%*******************************************
while 0
figure
FIG
plot(Cnt1(:,1), Cnt1(:,2))
plot(Cnt2(:,1), Cnt2(:,2), 'r')
end  % while 0
%*******************************************

nb_pt1 = size(Cnt1,1);
nb_pt2 = size(Cnt2,1);

% Double loop on all segment pairs
Pt_inter = [];
for i_pt1 = 1:nb_pt1-1
    for i_pt2 = 1:nb_pt2-1
        [X_inter, Y_inter] = inter_segs(Cnt1(i_pt1,1), Cnt1(i_pt1,2), Cnt1(i_pt1+1,1), Cnt1(i_pt1+1,2),...
                                        Cnt2(i_pt2,1), Cnt2(i_pt2,2), Cnt2(i_pt2+1,1), Cnt2(i_pt2+1,2));
        if ~isnan(X_inter)
            Pt_inter = [Pt_inter; X_inter, Y_inter];
        end  % if ~isnan(X_inter)
    end  % for i_pt2 = 1:nb_pt2-1
end  % for i_pt1 = 1:nb_pt1-1

%*******************************************
while 0
figure
FIG
plot(Cnt1(:,1), Cnt1(:,2))
plot(Cnt2(:,1), Cnt2(:,2), 'r')
plot(Pt_inter(:,1), Pt_inter(:,2), 'k*')
end  % while 0
%*******************************************


end












