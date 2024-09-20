% defines the RHS for the ode solver

function dx=RHS_fourvar(~,u,p)

    C = u(1:p.Nx);
    F = u(p.Nx+1:2*p.Nx);
    M = u(2*p.Nx+1:3*p.Nx);
    P = u(3*p.Nx+1:4*p.Nx);

    ghostleft = 2;
    ghostright = p.Nx - 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % n(i+1)
    xdiffR = [2:p.Nx, ghostright];
    % n(i-1)
    xdiffL = [ghostleft, 1:p.Nx-1];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End setting up vectors for differencing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Initialize the Jacobian as a sparse matrix
    Jac_block = sparse([1:p.Nx], [1:p.Nx], -2*ones(1, p.Nx), p.Nx, p.Nx);
    Jac_block = Jac_block + sparse([1:size(xdiffL,2)], xdiffL, ones(1, size(xdiffL,2)), p.Nx, p.Nx);
    Diffusion = Jac_block + sparse([1:size(xdiffR,2)], xdiffR, ones(1, size(xdiffR,2)), p.Nx, p.Nx);

%     Jac = p.Jac;
%     d = logical(eye(size(Jac)));
%     Jac(d) = Jac(d).*(-2);
%     
% 
%     C1 = p.Jac(p.Nx+1:2*p.Nx,1:p.Nx);
%     d1 = logical(eye(size(C1)));
%     C1(d1) = C1(d1).*(-2);
% 
%     C2 = p.Jac(p.Nx+1:2*p.Nx,1:p.Nx);
%     d2 = logical(eye(size(C2)));
%     C2(d2) = C2(d2).*0; 
%     [m, n] = size(C2);
%     lower_diag = tril(true(m, n), -1);
%     C2(lower_diag) = C2(lower_diag).*(-1);    
    
    diagval = ones(p.Nx,1);
    C1 = -2*spdiags(diagval, 0, p.Nx, p.Nx) + spdiags(diagval, -1, p.Nx, p.Nx) + spdiags(diagval, 1, p.Nx, p.Nx);
    C2 = -spdiags(diagval, -1, p.Nx, p.Nx) + spdiags(diagval, 1, p.Nx, p.Nx);

    Chemotaxis = (1/4).*((C2*F).*(C2*(1-C))) + (Diffusion*(1-C)).*F;


    dC = (1/(p.dx^2)).*Diffusion*C + p.sigma_C.*F.*(1-C).*(1+p.f_max.*M.*exp(1-M))-p.eta_bar.*M.*C-p.delta_C.*C;
    dF = (p.D_F/(p.dx^2)).*Diffusion*F - (p.chi_F/(p.dx^2)).*Chemotaxis + p.sigma_F.*F.*(1-F).*(1+(p.k_f.*C./(p.gamma_F+C))) - p.delta_F.*F;
    dM = (p.D_M/(p.dx^2)).*Diffusion*M + (p.k_M.*F).*(1-C)./(p.gamma_M+F) - p.mu_bar.*M.*P - p.delta_M.*M;
    dP = (p.D_P/(p.dx^2)).*Diffusion*P + p.k_P.*F.*M./(p.gamma_P+F) - p.mu_bar.*M.*P - p.delta_P.*P;
    
    dx	= [dC;dF;dM;dP];
