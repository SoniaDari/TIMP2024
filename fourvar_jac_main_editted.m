%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Simulation of four-variable model

    %% Problem set-up
    p.Nx = 1000;
    p.L = 100;
    x = linspace(0,p.L,p.Nx);
    p.dx = x(2)-x(1);
    q = 0;
    r = p.L;


    % Initial conditions
    C_init = heaviside(3-x);
    F_init = heaviside(3-x);
    M_init = 0.1*exp(-((x-3).^2));
    P_init = 0*ones(1,p.Nx);
    u0 = [C_init F_init M_init P_init];

    % set parameters
    D_f = 1.96*(10^(-6)); % Lit
    D_m = 9*(10^(-6)); % Lit
    D_p = 9*(10^(-6)); % Lit
    chi_f = 1.96*(10^(-6))*3;
    delta_c = 0.00001; % Estimate
    delta_f = 0.000139; % Lit
    delta_m = log(2)/49; % Lit
    delta_p = log(2)/4; % Lit
    
    sigma_c = 0.0183; % Lit
    c_0 = 0.68; % Lit
    p.f_max = 15; % Estimate
    sigma_f = 0.15;% Estimate
    f_0 = 0.632; % Lit
    gamma_m = 0.3; % Estimate
    k_p = 0.235; % % Lit
    gamma_p = 0.3; % Estimate

    eta = 9.144*10^(6); % Lit
    mu = 3.6*(10)^7; % Lit
    k_m = 2.523*10^(-10); % Lit
    m_thresh = 9*10^(-10); % Estimate
    
    p.sigma_C = f_0/c_0;
    p.eta_bar = eta*m_thresh/sigma_c;
    p.delta_C = delta_c/sigma_c;
%     p.delta_C = 2.5;
    p.chi_F = chi_f/D_f;
    p.sigma_F = sigma_f/sigma_c;
    p.delta_F = delta_f/sigma_c;
    p.D_M = D_m/D_f;
    p.k_M = k_m/(m_thresh*sigma_c);
%     p.k_M = 90;
    p.gamma_M = gamma_m/f_0;
    p.mu_bar = mu*m_thresh/sigma_c;
    p.delta_M = delta_m/sigma_c;
    p.D_P = D_p/D_f;
    p.k_P = k_p/sigma_c;
    p.gamma_P = gamma_p/f_0;
    p.delta_P = delta_p/sigma_c;
    
    tspan   = linspace(0,10,11);
    p.dt  = tspan(2)-tspan(1);
    nx = p.Nx;

    p.diffusion_coefficients = [0,1,p.D_M,p.D_P];
%% Define diffusion matrix & sparsity pattern for ode solver
    p.Jac = Jacobian_fourvar_editted(p);
%% ====================================================================
    % Run 1D code
options = odeset('Events',@wavespeed_events_fourvar);

[t,n,te,ye,ie]= ode15s(@RHS_fourvar_editted_reduced,tspan,u0,options,p);

    IE = ie<=nx;
    index = find(IE==1);
    te_ECM = te(index);
    if size(index) < 10
        e1 = 0;
    elseif size(index)<nx
        c1 = ((r-q)/nx)./diff(te_ECM);
        N = 3;
        NN = floor(length(c1)/N);
        d1 = c1(NN:2*NN);
        e1 = mean(d1);
    elseif size(index)<2*nx
        c1 = ((r-q)/nx)./diff(te_ECM(1:2:end));
        N = 3;
        NN = floor(length(c1)/N);
        d1 = c1(NN:2*NN);
        e1 = mean(d1);
    elseif size(index)<3*nx
        c1 = ((r-q)/nx)./diff(te_ECM(1:3:end));
        N = 3;
        NN = floor(length(c1)/N);
        d1 = c1(NN:2*NN);
        e1 = mean(d1);        
    else 
        e1 = 1;
    end
    sensitivity_metric = e1;


C = n(:,1:p.Nx);
F = n(:,p.Nx+1:2*p.Nx);
M = n(:,2*p.Nx+1:3*p.Nx);
P = n(:,3*p.Nx+1:4*p.Nx);

%% Figures
% figure(9)
% plot(x,C(1,:),'linewidth',2,'color','b')
% hold on
% for i = 1:11
%     hold on
%     plot(x,C(1*i,:),'linewidth',2,'color','b')
% end
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold')
% ylabel('C(X,T)','FontSize',18,'FontWeight','bold')
% ylim([0 1.2])
% xlim([0 100])
% 
% 
% figure(10)
% plot(x,F(1,:),'linewidth',2,'color','b')
% hold on
% for i = 1:10
%     hold on
%     plot(x,F(1*i,:),'linewidth',2,'color','b')
% end
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold')
% ylabel('F(X,T)','FontSize',18,'FontWeight','bold')
% ylim([0 1.2])
% xlim([0 100])
% 
% 
% figure(11)
% plot(x,M(1,:),'linewidth',2,'color','b')
% hold on
% for i = 1:10
%     hold on
%     plot(x,M(1*i,:),'linewidth',2,'color','b')
% end
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold')
% ylabel('M(X,T)','FontSize',18,'FontWeight','bold')
% xlim([0 100])
% 
% 
% figure(12)
% plot(x,P(1,:),'linewidth',2,'color','b')
% hold on
% for i = 1:10
%     hold on
%     plot(x,P(1*i,:),'linewidth',2,'color','b')
% end
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold')
% ylabel('P(X,T)','FontSize',18,'FontWeight','bold')
% xlim([0 100])




    