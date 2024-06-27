% [Tupdt] = repl_val(T, val1, val2);
% [T1updt, T2updt] = repl_val(T1, T2, val1, val2);
% [T1updt, T2updt, T3updt] = repl_val(T1, T2, T3, val1, val2);
% 
% Replaces in T all the values val1 by val2.
% If several Ts are considered, replacement occurs only if the value to
% replace exists in the multiple Ts
% 
% Inputs
%   T, T1, T2, T3 : Input matrices (of same size)
%   val1          : Value to replace
%   val2          : Replacement value
% 
% Outputs
%   T, T1, T2, T3 : Output matrices
% 
% Author (adaptation): Antoine Serrurier
% Date (adaptation): 25/06/2024

function [T1_cor, T2_cor, T3_cor] = repl_val(P1, P2, P3, P4, P5);

% Check the number of orguments
if ~(nargin == 3 & nargout == 1) & ~(nargin == 4 & nargout == 2) & ~(nargin == 5 & nargout == 3)
  error('Error in repl_val: wrong number of input arguments')
end
  
% Dimensions and renaming input arguments
if nargin == 3
  T1 = P1; val1 = P2; val2 = P3;
  [T1_n1, T1_n2] = size(T1);
elseif nargin == 4
  T1 = P1; T2 = P2; val1 = P3; val2 = P4;
  [T1_n1, T1_n2] = size(T1);
  [T2_n1, T2_n2] = size(T2);
elseif nargin == 5
  T1 = P1; T2 = P2; T3 = P3; val1 = P4; val2 = P5;
  [T1_n1, T1_n2] = size(T1);
  [T2_n1, T2_n2] = size(T2);
  [T3_n1, T3_n2] = size(T3);
end

% Check val1
if size(val1, 1)~=1 | size(val1, 2)~=1
  error('Error in repl_val: val1 must be a simple value')
end

if size(val2, 1)~=1 | size(val2, 2)~=1
  error('Error in repl_val: val2 must be a simple value')
end

% Treat each case one after the other

if nargin == 3 % -------------------------------------------------
T1_cor = T1;
if isnan(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isnan(T1(i, j))
	T1_cor(i, j) = val2;
      end
    end
  end

elseif isinf(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isinf(T1(i, j))
	T1_cor(i, j) = val2;
      end
    end
  end

else
  for i = 1:T1_n1
    for j = 1:T1_n2
      if T1(i, j) == val1
	T1_cor(i, j) = val2;
      end
    end
  end
end

elseif nargin == 4 % -------------------------------------------------
T1_cor = T1;
T2_cor = T2;
if isnan(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isnan(T1(i, j)) & isnan(T2(i, j))
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
      end
    end
  end

elseif isinf(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isinf(T1(i, j)) & isinf(T2(i, j))
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
      end
    end
  end

else
  for i = 1:T1_n1
    for j = 1:T1_n2
      if T1(i, j) == val1 & T2(i, j) == val1
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
      end
    end
  end
end

elseif nargin == 5 % -------------------------------------------------
T1_cor = T1;
T2_cor = T2;
T3_cor = T3;
if isnan(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isnan(T1(i, j)) & isnan(T2(i, j)) & isnan(T3(i, j))
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
	T3_cor(i, j) = val2;
      end
    end
  end

elseif isinf(val1)
  for i = 1:T1_n1
    for j = 1:T1_n2
      if isinf(T1(i, j)) & isinf(T2(i, j)) & isinf(T3(i, j))
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
	T3_cor(i, j) = val2;
      end
    end
  end

else
  for i = 1:T1_n1
    for j = 1:T1_n2
      if T1(i, j) == val1 & T2(i, j) == val1 & T3(i, j) == val1
	T1_cor(i, j) = val2;
	T2_cor(i, j) = val2;
	T3_cor(i, j) = val2;
      end
    end
  end
end

end % elseif nargin == 5 % -------------------------------------------------

return
