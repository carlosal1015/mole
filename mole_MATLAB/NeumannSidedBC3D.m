function BC = NeumannSidedBC3D(k, m, dx, left, left_coef, ...
                                        right, right_coef, ...
                                        front, front_coef, ...
                                         back, back_coef, ...
                                          top, top_coef, ...
                                       bottom, bottom_coef)
% Returns a three-dimensional mimetic boundary operator that 
% imposes a boundary condition of Robin's type
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells along x-axis
%               dx : Step size along x-axis
%                n : Number of cells along y-axis
%               dy : Step size along y-axis
%                o : Number of cells along z-axis
%               dz : Step size along z-axis
%              left: T/F (binary) build left BC
%        left_coef : Value for left BC
%             right: T/F (binary) build right BC
%       right_coef : Value for right BC
%             front: T/F (binary) build front BC
%       front_coef : Value for front BC
%             back : T/F (binary) build back BC
%        back_coef : Value for back BC
%              top : T/F (binary) build top BC
%         top_coef : Value for top BC
%           bottom : T/F (binary) build bottom BC
%      bottom_coef : Value for bottom BC


    % 1-D boundary operator
    Bm = NeumannSidedBC(k, m, dx, left, left_coef, right, right_coef);
    Bn = NeumannSidedBC(k, n, dy, front, front_coef, back, back_coef);
    Bo = NeumannSidedBC(k, o, dz, top, top_coef, bottom, bottom_coef);
    
    Im = speye(m+2);
    In = speye(n+2);
    Io = speye(o+2);
    
    Io(1, 1) = 0;
    Io(end, end) = 0;
    
    In2 = In;
    In2(1, 1) = 0;
    In2(end, end) = 0;
    
    BC1 = kron(kron(Io, In2), Bm);
    BC2 = kron(kron(Io, Bn), Im);
    BC3 = kron(kron(Bo, In), Im);
    
    BC = BC1 + BC2 + BC3;
end
