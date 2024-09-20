function Jac=Jacobian_fourvar(p)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  sets up difference matrix
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  no-flux boundary conditions
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ghostleft  = 2;
% ghostright = p.Nx-1; 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  n(i+1)
% xdiffR	   = [2:p.Nx,ghostright];
% %  n(i-1)
% xdiffL	   = [ghostleft,1:p.Nx-1];
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % End setting up vectors for differencing 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %  Diffusion terms:
% % D(n_i+/-1 - 2 n_i)
% Jac=sparse([1:p.Nx],[1:p.Nx],-2*ones(1,p.Nx),p.Nx,p.Nx);
% %  
% Jac=Jac+sparse([1:size(xdiffL,2)],xdiffL,ones(1,size(xdiffL,2)),p.Nx,p.Nx);
% Jac=Jac+sparse([1:size(xdiffR,2)],xdiffR,ones(1,size(xdiffR,2)),p.Nx,p.Nx);
% 
% 
% % Diffusion terms for a single equation
% baseBlock = sparse([1:p.Nx], [1:p.Nx], -2 * ones(1, p.Nx), p.Nx, p.Nx);
% Jac = kron(speye(4), baseBlock);
% 
% % Update the diagonal blocks with diffusion terms
% Jac(sub2ind(size(Jac), kron(1:4, ones(1, p.Nx)), xdiffL)) = 1;
% Jac(sub2ind(size(Jac), kron(1:4, ones(1, p.Nx)), xdiffR)) = 1;
% 
% clear xdiffL xdiffR ghostleft ghostright
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sets up difference matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no-flux boundary conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     ghostleft = 2;
%     ghostright = p.Nx - 1;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % n(i+1)
%     xdiffR = [2:p.Nx, ghostright];
%     % n(i-1)
%     xdiffL = [ghostleft, 1:p.Nx-1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End setting up vectors for differencing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialize the Jacobian as a sparse matrix
%     Jac_block = sparse([1:p.Nx], [1:p.Nx], ones(1, p.Nx), p.Nx, p.Nx);
%     Jac_block = Jac_block + sparse([1:size(xdiffL,2)], xdiffL, ones(1, size(xdiffL,2)), p.Nx, p.Nx);
%     Jac_block = Jac_block + sparse([1:size(xdiffR,2)], xdiffR, ones(1, size(xdiffR,2)), p.Nx, p.Nx);

    fill = ones(p.Nx,1);
    fill_mat = spdiags([fill,fill,fill],-1:1,p.Nx,p.Nx);

    
% %     D_C = Jac_block./p.dx^2;
% %     D_F = p.D_F.*Jac_block./p.dx^2;
% %     D_M = p.D_M.*Jac_block./p.dx^2;
% %     D_P = p.D_P.*Jac_block./p.dx^2;
%     Jac = blkdiag(Jac_block,D_F,D_M,D_P);
%     spy(Jac);

% % Chemotaxis
%     C_middle = sparse([1:p.Nx], [1:p.Nx], ones(1, p.Nx), p.Nx, p.Nx);
%     C_lower = sparse([1:size(xdiffL,2)], xdiffL, ones(1, size(xdiffL,2)), p.Nx, p.Nx);
%     C_upper = sparse([1:size(xdiffR,2)], xdiffR, ones(1, size(xdiffR,2)), p.Nx, p.Nx);
% 
%     F_middle = Jac_block;
%     F_lower = -sparse([1:p.Nx], [1:p.Nx], ones(1, p.Nx), p.Nx, p.Nx) + sparse([1:size(xdiffL,2)], xdiffL, ones(1, size(xdiffL,2)), p.Nx, p.Nx);
%     F_upper = -sparse([1:p.Nx], [1:p.Nx], ones(1, p.Nx), p.Nx, p.Nx) + sparse([1:size(xdiffR,2)], xdiffR, ones(1, size(xdiffR,2)), p.Nx, p.Nx);
% 
%     Chemotaxis = C_middle*F_middle + C_lower*F_lower + C_upper*F_upper;
%     Chemotaxis = Chemotaxis*p.chi_F/(2*(p.dx^2));


    id = (sparse([1:p.Nx], [1:p.Nx], ones(1, p.Nx), p.Nx, p.Nx));
%     id = logical(ones(p.Nx,p.Nx));
    zero_mat = id;

    
    row1 = [fill_mat,zero_mat,zero_mat,zero_mat];
    row2 = [fill_mat,fill_mat,zero_mat,zero_mat];
    row3 = [zero_mat,zero_mat,fill_mat,zero_mat];
    row4 = [zero_mat,zero_mat,zero_mat,fill_mat];

    Jac = [row1;row2;row3;row4];    
    
    
    % Clear temporary variables
    clear xdiffL xdiffR ghostleft ghostright
end






