%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Simulation of four-variable model

    %% Problem set-up
    p.Nx = 100;
    p.L = 100;
    x = linspace(0,p.L,p.Nx);
    p.dx = x(2)-x(1);
    % Initial conditions
    C_init = heaviside(3-x);
    F_init = heaviside(3-x);
    M_init = 0.1*exp(-((x-3).^2));
    P_init = 0*ones(1,p.Nx);
    u0 = [C_init F_init M_init P_init];

    % set parameters
    D_c = 0.05882352941*6.12*(10^(-7));
    D_f = 6.12*(10^(-7));
    D_m = 9*(10^(-6)); 
    D_p = 9*(10^(-6))/37;
    chi_f = 6.12*(10^(-7))*10/3;
%     chi_f = 0;
    delta_c = 0.001; % Estimate
    delta_f = 0.005; 
    delta_m = 0.009; 
    delta_p = 0.1925;
    
    sigma_c = 0.018;
    c_0 = 0.1; % Estimate
    p.f_max = 50; % Estimate
    p.k_f = 30; %Estimate
    gamma_f = 0.05; % Estimate
    m_thresh = 5*10^(-9); % Estimate
    eta = 1.08*10^(6);
    sigma_f = 0.0385;
    f_0 = 0.12;
    k_m = 4*2.9*10^(-10);
    gamma_m = 0.05; % Estimate
    mu = 4.98*10^(8)/24; 
    k_p = 3.86*10^(-2); % Estinate
    gamma_p = 0.05; % Estimate
    
    p.sigma_C = f_0/c_0;
    p.eta_bar = eta*m_thresh/sigma_c;
    p.delta_C = delta_c/sigma_c;
    p.D_F = D_f/D_c;
    p.chi_F = chi_f/D_c;
    p.sigma_F = sigma_f/sigma_c;
    p.gamma_F = gamma_f/c_0;
    p.delta_F = delta_f/sigma_c;
    p.D_M = D_m/D_c;
    p.k_M = k_m/(m_thresh*sigma_c);
    p.gamma_M = gamma_m/f_0;
    p.mu_bar = mu*m_thresh/sigma_c;
    p.delta_M = delta_m/sigma_c;
    p.D_P = D_p/D_c;
    p.k_P = k_p/sigma_c;
    p.k_P = 0;
    p.gamma_P = gamma_p/f_0;
    p.delta_P = delta_p/sigma_c;

    p.diffusion_coefficients = [1,p.D_F,p.D_M,p.D_P];
%% Define diffusion matrix & sparsity pattern for ode solver
%     p.Jac = Jacobian_fourvar(p)./p.dx^2 ;
    p.Jac = Jacobian_fourvar(p);
%% ====================================================================
    % Run 1D code
%     options     = odeset('AbsTol',1e-5,'RelTol',1e-5,'JPattern',p.Jac,'Events',@wave_events);
%     [t,y,TE,~,IE]= ode15s(@pushpull_equations,[0 tend],xzero,options,p);% time points as wave crosses each node
%     c           = 1./diff(TE(2:end));                                   % calculate speed

%%
tspan   = linspace(0,20,20);
p.dt  = tspan(2)-tspan(1);
% options = odeset('JPattern',p.Jac);
% [t u] 	= ode15s(@RHS_fourvar,tspan,u0,options,p);

    options     = odeset('AbsTol',1e-5,'RelTol',1e-5,'JPattern',p.Jac,'Events',@wavespeed_events_fourvar);
    [t,n,te,ye,ie]= ode15s(@RHS_fourvar,tspan,u0,options,p);% time points as wave crosses each node
%     c           = 1./diff(TE(2:end));                                   % calculate speed


C = n(:,1:p.Nx);
F = n(:,p.Nx+1:2*p.Nx);
M = n(:,2*p.Nx+1:3*p.Nx);
P = n(:,3*p.Nx+1:4*p.Nx);

%% Figures
figure(5)
plot(x,C(1,:),'linewidth',2,'color','k')
hold on
for i = 1:20
    hold on
    plot(x,C(1*i,:),'linewidth',2,'color','k')
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('C(X,T)','FontSize',18,'FontWeight','bold')


figure(6)
plot(x,F(1,:),'linewidth',2,'color','k')
hold on
for i = 1:20
    hold on
    plot(x,F(1*i,:),'linewidth',2,'color','k')
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('F(X,T)','FontSize',18,'FontWeight','bold')


figure(7)
plot(x,M(1,:),'linewidth',2,'color','k')
hold on
for i = 1:20
    hold on
    plot(x,M(1*i,:),'linewidth',2,'color','k')
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('M(X,T)','FontSize',18,'FontWeight','bold')


figure(8)
plot(x,P(1,:),'linewidth',2,'color','k')
hold on
for i = 1:20
    hold on
    plot(x,P(1*i,:),'linewidth',2,'color','k')
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold')
ylabel('P(X,T)','FontSize',18,'FontWeight','bold')




    