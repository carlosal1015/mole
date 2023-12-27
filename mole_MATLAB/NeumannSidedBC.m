function BC = NeumannSidedBC(k, m, dx, left, left_coef, right, right_coef)
% Returns a m+2 by m+2 one-dimensional mimetic boundary operator that 
% imposes a boundary condition of Neumann  type
%
% Parameters:
%                k : Order of accuracy
%                m : Number of cells
%               dx : Step size
%              left: T/F (binary) build left BC
%        left_coef : Value for left BC
%             right: T/F (binary) build right BC
%       right_coef : Value for right BC
    
    % Dirichlet could be made in a similar fashion
    %A = sparse(m+2, m+2);
    
    B = sparse(m+2, m+1);
    
    if ( left )
        B(1, 1) = -left_coef;
    end
    
    if ( right )
        B(end, end) = right_coef;
    end
    
    G = grad(k, m, dx);
    
    BC = B*G;
    
    % Or if making Robin/Dirichlet
    %BC = A + B*G;
end
