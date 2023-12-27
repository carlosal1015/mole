function BC = NeumannSidedBC2D(k, m, dx, n, dy, left, left_coef, right, right_coef, front, front_coef, back, back_coef)
% Returns a two-dimensional mimetic boundary operator that 
% imposes a boundary condition of Robin's type
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%              left: T/F (binary) build left BC
%        left_coef : Value for left BC
%            right : T/F (binary) build right BC
%       right_coef : Value for right BC
%            front : T/F (binary) build front BC
%       front_coef : Value for front(bottom) BC
%              back: T/F (binary) build back BC
%        back_coef : Value for back(top) BC

    % 1-D boundary operator
    Bm = NeumannSidedBC(k, m, dx, left, left_coef, right, right_coef);
    Bn = NeumannSidedBC(k, n, dy, front, front_coef, back, back_coef);
    
    Im = speye(m+2);
    In = speye(n+2);
    
    In(1, 1) = 0;
    In(end, end) = 0;
    
    BC1 = kron(In, Bm);
    BC2 = kron(Bn, Im);
    
    BC = BC1 + BC2;
end
