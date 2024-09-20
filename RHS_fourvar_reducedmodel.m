% defines the RHS for the ode solver

function dx=RHS_fourvar_reducedmodel(~,u,p)

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

    dC = (1/(p.dx^2)).*Diffusion*C + p.sigma_C.*F.*(1-C).*(1+p.f_max.*M.*exp(1-M))-p.eta_bar.*M.*C-p.delta_C.*C;
    dF = (p.D_F/(p.dx^2)).*Diffusion*F + p.sigma_F.*F.*(1-F).*(1+(p.k_f.*C)) - p.delta_F.*F;
    dM = (p.D_M/(p.dx^2)).*Diffusion*M + (p.k_M.*F).*(1-C) - p.mu_bar.*M.*P - p.delta_M.*M;
    dP = (p.D_P/(p.dx^2)).*Diffusion*P + p.k_P.*F.*M - p.mu_bar.*M.*P - p.delta_P.*P;
    
    dx	= [dC;dF;dM;dP];
