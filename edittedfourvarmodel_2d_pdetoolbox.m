%% Four variable model with nonlinearity
%% Parameter Values
% Dimensional
    D_f = 1.96*(10^(-6)); % Lit
    D_m = 9*(10^(-6)); % Lit
    D_p = 9*(10^(-6)); % Lit
    delta_c = 0.00001; % Estimate
    delta_f = 0.000139; % Lit
    delta_m = log(2)/49; % Lit
    delta_p = log(2)/4; % Lit
    
    sigma_c = 0.0183; % Lit
    c_0 = 0.68; % Lit
    f_max = 15; % Estimate
    sigma_f = 0.15;% Estimate
    f_0 = 0.632; % Lit
    k_p = 0.235; % % Lit

    eta = 9.144*10^(6); % Lit
    mu = 3.6*(10)^7; % Lit
    k_m = 2.523*10^(-10); % Lit
    m_thresh = 9*10^(-10); % Estimate
    
    sigma_C = f_0/c_0;
    eta_bar = eta*m_thresh/sigma_c;
    delta_C = delta_c/sigma_c;
    sigma_F = sigma_f/sigma_c;
    delta_F = delta_f/sigma_c;
    D_M = D_m/D_f;
    k_M = k_m/(m_thresh*sigma_c);
    mu_bar = mu*m_thresh/sigma_c;
    delta_M = delta_m/sigma_c;
    D_P = D_p/D_f;
    k_P = k_p/sigma_c;
    delta_P = delta_p/sigma_c;

%% Initialise PDE
model = createpde(4);
% model.SolverOptions.RelativeTolerance = 1E-7;
% model.SolverOptions.AbsoluteTolerance = 1E-11;
L=5;
M=100;
%% Polygonal Domain
poly1 = [3
    4
    0
    0
    M
    M
    0
    L
    L
    0
    ];

geom = [poly1];
ns = char('poly1');
ns = ns';
sf = 'poly1';
gd = decsg(geom,sf,ns);

%% Plot Computational Domain
fig = figure(1)
pdegplot(gd,'Edgelabels','on','Facelabels','on')
axis equal
% saveas(fig,'diffdomain.eps','epsc')

%% Attach domain to PDE model
geometryFromEdges(model,gd);

%% Set boundary conditions
% Zero flux BCs
applyBoundaryCondition(model,'neumann','edge',[1,2,3,4],'g',0,'q',0);

% Robin BC to model influx of TIMPs from hydrogel on side 2 (top)
% q2 = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1];
% g2 = [0 0 0 100];
% applyBoundaryCondition(model,'neumann','edge',[1,3,4],'q',0,'g',0);
% applyBoundaryCondition(model,'neumann','edge',2,'q',q2,'g',g2);


%% Specify coefficients for PDE model
m = 0;
d = [1;1;1;1];
a = [0;0;0;0];
c = [0.01;0.01;1;1;D_M;D_M;D_P;D_P];

specifyCoefficients(model,'m',m,'d',d,'c',c,'a',a,'f',@fcoeffun);

%% Generate the finite element mesh
m1=generateMesh(model,"Hmax",0.5);

%% Plot mesh
figure(2)
pdeplot(m1)
axis equal

%% Establish time domain
tlist = 0:1:20;

%% Set initial conditions
setInitialConditions(model,@ut0fun)

%% Solve the PDE
result = solvepde(model,tlist);

% Extract solution data
u = result.NodalSolution;



%% Visulaisation
fig = figure(1);
% Initialise VideoWriter object to capture a movie of the simulations
vid0bj4varmodel_C = VideoWriter('fourvarmodelnonlinearity2D_C.avi');
open(vid0bj4varmodel_C);
for i = 1:length(tlist)
    pdeplot(model,'XYData',u(:,1,i),'ZData',u(:,1,i),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
    title(strjoin(["C(X,Y,T); T=",num2str(tlist(i))]));
    xlabel('X')
    ylabel('Y')
    zlabel('C(X,Y,T)')
    zlim([0,1.1])
    drawnow limitrate
    F(i) = getframe(fig);
    writeVideo(vid0bj4varmodel_C,F(i));
end
close(vid0bj4varmodel_C);

fig2 = figure(2);
% Initialise VideoWriter object to capture a movie of the simulations
vid0bj4varmodel_F = VideoWriter('fourvarmodelnonlinearity2D_F.avi');
open(vid0bj4varmodel_F);
for i = 1:length(tlist)
    pdeplot(model,'XYData',u(:,2,i),'ZData',u(:,2,i),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
    title(strjoin(["F(X,Y,T); T=",num2str(tlist(i))]));
    xlabel('X')
    ylabel('Y')
    zlabel('F(X,Y,T)')
    zlim([0,1.1])
    drawnow limitrate
    F(i) = getframe(fig2);
    writeVideo(vid0bj4varmodel_F,F(i));
end
close(vid0bj4varmodel_F);

fig3 = figure(3);
% Initialise VideoWriter object to capture a movie of the simulations
vid0bj4varmodel_M = VideoWriter('fourvarmodelnonlinearity2D_M.avi');
open(vid0bj4varmodel_M);
for i = 1:length(tlist)
    pdeplot(model,'XYData',u(:,3,i),'ZData',u(:,3,i),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
    title(strjoin(["M(X,Y,T); T=",num2str(tlist(i))]));
    xlabel('X')
    ylabel('Y')
    zlabel('M(X,Y,T)')
    drawnow limitrate
    F(i) = getframe(fig3);
    writeVideo(vid0bj4varmodel_M,F(i));
end
close(vid0bj4varmodel_F);

fig4 = figure(4);
% Initialise VideoWriter object to capture a movie of the simulations
vid0bj4varmodel_P = VideoWriter('fourvarmodelnonlinearity2D_P.avi');
open(vid0bj4varmodel_P);
for i = 1:length(tlist)
    pdeplot(model,'XYData',u(:,4,i),'ZData',u(:,4,i),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
    title(strjoin(["P(X,Y,T); T=",num2str(tlist(i))]));
    xlabel('X')
    ylabel('Y')
    zlabel('P(X,Y,T)')
    drawnow limitrate
    F(i) = getframe(fig4);
    writeVideo(vid0bj4varmodel_P,F(i));
end
close(vid0bj4varmodel_P);

% xq = linspace(0,100,10000);
% yq = ones(size(xq));
% component = 3;
% uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
% uintrp = squeeze(uintrp);

% [X,Y] = ndgrid(xq,tlist);
% figure(10)
% surf(X,Y,uintrp)
% xlabel("x")
% ylabel("Time")
% hold on

%%
% figure(5)
% pdeplot(model,'XYData',u(:,1,11),'ZData',u(:,1,11),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
% ylabel('Y','FontSize',18,'FontWeight','bold','interpreter','latex')
% zlabel('$C(X,Y,T_0)$','FontSize',18,'FontWeight','bold','interpreter','latex')
% xlim([0 100])
% box on
% grid on
% 
% figure(6)
% pdeplot(model,'XYData',u(:,2,11),'ZData',u(:,2,11),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
% ylabel('Y','FontSize',18,'FontWeight','bold','interpreter','latex')
% zlabel('$F(X,Y,T_0)$','FontSize',18,'FontWeight','bold','interpreter','latex')
% xlim([0 100])
% box on
% grid on
% 
% figure(7)
% pdeplot(model,'XYData',u(:,3,11),'ZData',u(:,3,11),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
% ylabel('Y','FontSize',18,'FontWeight','bold','interpreter','latex')
% zlabel('$M(X,Y,T_0)$','FontSize',18,'FontWeight','bold','interpreter','latex')
% xlim([0 100])
% box on
% grid on
% 
% figure(8)
% pdeplot(model,'XYData',u(:,4,11),'ZData',u(:,4,11),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% set(gca,'fontsize',16)
% xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
% ylabel('Y','FontSize',18,'FontWeight','bold','interpreter','latex')
% zlabel('$P(X,Y,T_0)$','FontSize',18,'FontWeight','bold','interpreter','latex')
% xlim([0 100])
% box on
% grid on

%%
% xq = linspace(0,100,1000);
% yq_1 = 1*ones(size(xq));
% yq = yq_1;
% component = 4;
% uintrp = interpolateSolution(result,xq,yq_1,component,1:length(tlist));
% uintrp_1 = squeeze(uintrp);
% 
% xq = linspace(0,100,1000);
% yq_2 = 2*ones(size(xq));
% component = 4;
% uintrp = interpolateSolution(result,xq,yq_2,component,1:length(tlist));
% uintrp_2 = squeeze(uintrp);
% 
% xq = linspace(0,100,1000);
% yq_3 = 3*ones(size(xq));
% component = 4;
% uintrp = interpolateSolution(result,xq,yq_3,component,1:length(tlist));
% uintrp_3 = squeeze(uintrp);
% 
% xq = linspace(0,100,1000);
% yq_4 = 4*ones(size(xq));
% component = 4;
% uintrp = interpolateSolution(result,xq,yq_4,component,1:length(tlist));
% uintrp_4 = squeeze(uintrp);
% 
% xq = linspace(0,100,1000);
% yq_5 = 5*ones(size(xq));
% component = 4;
% uintrp = interpolateSolution(result,xq,yq_5,component,1:length(tlist));
% uintrp_5 = squeeze(uintrp);
% 
% matrices = cat(3, uintrp_1, uintrp_2, uintrp_3, uintrp_4, uintrp_5);
% avgMatrix = mean(matrices, 3);
% 
% 
% Num_x = length(xq);
% N_avg = 1000;
% yq = ones(N_avg,Num_x);
% rows = linspace(0,5,N_avg);
% uintrp = zeros(Num_x,length(tlist),N_avg);
% 
% for i=1:N_avg
%     yq(i,:) = yq(i,:).*rows(i);
% end
% 
% for i=1:N_avg
%     uintrp(:,:,i) = interpolateSolution(result,xq,yq(i,:),component,1:length(tlist));
% end
% 
% avgMatrix = mean(unitrp, 3);
% 
% 
% 
% 
% 

% xq = linspace(0,100,1001);
% yq = ones(size(xq));
% component = 1;
% 
% uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
% uintrp = squeeze(uintrp);
% plot(xq,uintrp)


xq = linspace(0,100,10000);
yq = 0.05*ones(size(xq));
component = 1;
uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
uintrp = squeeze(uintrp);
figure(9)
plot(xq,heaviside(3-xq),'color','k','linewidth',2)
% plot(xq,heaviside(3-xq) + heaviside(xq-97),'color','k','linewidth',2)
for i = 1:10
hold on
plot(xq,uintrp(:,2*i),'color','k','linewidth',2)
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
ylabel('$C(X,T)$','FontSize',18,'FontWeight','bold','interpreter','latex')
box on
ylim([0,1.2])

component = 2;
uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
uintrp = squeeze(uintrp);
figure(10)
plot(xq,heaviside(3-xq),'color','k','linewidth',2)
plot(xq,heaviside(3-xq) + heaviside(xq-97),'color','k','linewidth',2)
for i = 1:10
hold on
plot(xq,uintrp(:,2*i),'color','k','linewidth',2)
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
ylabel('$F(X,\overline{Y},T)$','FontSize',18,'FontWeight','bold','interpreter','latex')
box on
ylim([0,1.2])

component = 3;
uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
uintrp = squeeze(uintrp);
figure(11)
plot(xq,0.1*exp(-((xq-3).^2)),'color','k','linewidth',2)
plot(xq,0.1*exp(-((xq-3).^2)) + 0.1*exp(-((xq-97).^2)),'color','k','linewidth',2)
for i = 1:10
hold on
plot(xq,uintrp(:,2*i),'color','k','linewidth',2)
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
ylabel('$M(X,\overline{Y},T)$','FontSize',18,'FontWeight','bold','interpreter','latex')
box on

component = 4;
uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
uintrp = squeeze(uintrp);
figure(12)
plot(xq,0*ones(length(xq),1),'color','k','linewidth',2)
for i = 1:10
hold on
plot(xq,uintrp(:,2*i),'color','k','linewidth',2)
end
set(gca,'fontsize',16)
xlabel('X','FontSize',18,'FontWeight','bold','interpreter','latex')
ylabel('$P(X,\overline{Y},T)$','FontSize',18,'FontWeight','bold','interpreter','latex')
box on

%%
% f = figure;
% f.Position = [100 100 1000 600];
% hold on
% pos1 = [0.2 0.75 0.13 0.17];
% subplot('Position',pos1)
% pdeplot(model,'XYData',u(:,1,9),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylabel('Y','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'fontsize',10);
% box on
% 
% pos2 = [0.35 0.75 0.13 0.17];
% subplot('Position',pos2)
% pdeplot(model,'XYData',u(:,1,17),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off') 
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos3 = [0.5 0.75 0.13 0.17];
% subplot('Position',pos3)
% pdeplot(model,'XYData',u(:,1,25),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos4 = [0.65 0.75 0.13 0.17];
% subplot('Position',pos4)
% pdeplot(model,'XYData',u(:,1,33),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos5 = [0.8 0.75 0.18 0.17];
% subplot('Position',pos5)
% pdeplot(model,'XYData',u(:,1,41),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% set(gca,'xtick',[],'ytick',[],'fontsize',10);
% box on
% 
% 
% 
% 
% pos6 = [0.2 0.55 0.13 0.17];
% subplot('Position',pos6)
% pdeplot(model,'XYData',u(:,2,9),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylabel('Y','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'fontsize',10);
% box on
% 
% pos7 = [0.35 0.55 0.13 0.17];
% subplot('Position',pos7)
% pdeplot(model,'XYData',u(:,2,17),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off') 
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos8 = [0.5 0.55 0.13 0.17];
% subplot('Position',pos8)
% pdeplot(model,'XYData',u(:,2,25),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos9 = [0.65 0.55 0.13 0.17];
% subplot('Position',pos9)
% pdeplot(model,'XYData',u(:,2,33),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos10 = [0.8 0.55 0.18 0.17];
% subplot('Position',pos10)
% pdeplot(model,'XYData',u(:,2,41),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 1])
% set(gca,'xtick',[],'ytick',[],'fontsize',10);
% box on
% 
% 
% 
% 
% 
% 
% pos11 = [0.2 0.35 0.13 0.17];
% subplot('Position',pos11)
% pdeplot(model,'XYData',u(:,3,9),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylabel('Y','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 12])
% colorbar('off')
% set(gca,'xtick',[],'fontsize',10);
% box on
% 
% pos12 = [0.35 0.35 0.13 0.17];
% subplot('Position',pos12)
% pdeplot(model,'XYData',u(:,3,17),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 12])
% colorbar('off') 
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos13 = [0.5 0.35 0.13 0.17];
% subplot('Position',pos13)
% pdeplot(model,'XYData',u(:,3,25),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 12])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos14 = [0.65 0.35 0.13 0.17];
% subplot('Position',pos14)
% pdeplot(model,'XYData',u(:,3,33),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 12])
% colorbar('off')
% set(gca,'xtick',[],'ytick',[]);
% box on
% 
% pos15 = [0.8 0.35 0.18 0.17];
% subplot('Position',pos15)
% pdeplot(model,'XYData',u(:,3,41),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylim([0 5])
% xlim([0 100])
% caxis([0 12])
% set(gca,'xtick',[],'ytick',[],'fontsize',10);
% box on
% 
% 
% 
% 
% 
% pos16 = [0.2 0.15 0.13 0.17];
% subplot('Position',pos16)
% pdeplot(model,'XYData',u(:,4,9),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% ylabel('Y','FontSize',12)
% xlabel('X','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 50])
% colorbar('off')
% set(gca,'fontsize',10);
% box on
% 
% pos17 = [0.35 0.15 0.13 0.17];
% subplot('Position',pos17)
% pdeplot(model,'XYData',u(:,4,17),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% xlabel('X','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 50])
% colorbar('off') 
% set(gca,'ytick',[],'fontsize',10);
% box on
% 
% pos18 = [0.5 0.15 0.13 0.17];
% subplot('Position',pos18)
% pdeplot(model,'XYData',u(:,4,25),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% xlabel('X','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 50])
% colorbar('off')
% set(gca,'ytick',[],'fontsize',10);
% box on
% 
% pos19 = [0.65 0.15 0.13 0.17];
% subplot('Position',pos19)
% pdeplot(model,'XYData',u(:,4,33),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% xlabel('X','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 50])
% colorbar('off')
% set(gca,'ytick',[],'fontsize',10);
% box on
% 
% pos20 = [0.8 0.15 0.18 0.17];
% subplot('Position',pos20)
% pdeplot(model,'XYData',u(:,4,41),'FaceAlpha',0.5,'ColorMap','parula','Mesh','off')
% xlabel('X','FontSize',12)
% ylim([0 5])
% xlim([0 100])
% caxis([0 50])
% set(gca,'ytick',[],'fontsize',10);
% box on





%% Wave Speed
% xq = linspace(0,100,10000);
% yq = ones(size(xq));
% x = xq;
% t = tlist;
% component = 1;
% uintrp = interpolateSolution(result,xq,yq,component,1:length(tlist));
% uintrp = squeeze(uintrp);
% U = uintrp;
% a = length(U);
% b = width(U);
% % Assuming U is your a x b matrix where a represents spatial points and b represents time points
% % Also assuming x is your spatial points vector and t is your time points vector
% 
% % Define the threshold value for identifying the wavefront
% threshold = 0.5; % Set an appropriate threshold value
% 
% % Initialize an array to store the wavefront positions
% % wavefront_positions = zeros(3, 16);
% 
% % Preallocate wavefront_positions
% wavefront_positions = zeros(1, 16); % Assuming 16 time points
% 
% % Find the wavefront position for each time point
% for j = 1:16
%     % Find the indices where the condition is met
%     A = find(U(:, j) <= threshold);
%     
%     if ~isempty(A)
%         % Get the first index where U crosses the threshold
%         wavefront_index = min(A);
%         wavefront_positions(j) = x(wavefront_index);
%     else
%         % If no value meets the condition, assign NaN or some placeholder
%         wavefront_positions(j) = NaN; % No crossing found
%     end
% end
% 
% 
% % Calculate the displacements of the wavefront over time
% displacements = diff(wavefront_positions);
% 
% % Calculate the time intervals
% time_intervals = diff(t);
% 
% % Calculate the wave speed for each time interval
% wave_speeds = displacements ./ time_intervals(1);
% 
% % Calculate the average wave speed
% average_wave_speed = mean(wave_speeds, 'omitnan');
% 
% % Display the result
% fprintf('The approximate wave speed is %.4f units/time\n', average_wave_speed);


%% Function Declarations
function ut0 = ut0fun(location)
    nr = length(location.x);
    ut0(1,:) = heaviside(3-location.x) + heaviside(location.x-97);
    ut0(2,:) = heaviside(3-location.x) + heaviside(location.x-97);
    ut0(3,:) = 0.1*exp(-((location.x-3).^2)) + 0.1*exp(-((location.x-97).^2));
    ut0(4,:) = 0;
end

% function ut0 = ut0fun(location)
%     nr = length(location.x);
%     ut0(1,:) = heaviside(3-location.x);
%     ut0(2,:) = heaviside(3-location.x);
%     ut0(3,:) = 0.1*exp(-((location.x-3).^2));
%     ut0(4,:) = 0;
% end


function f = fcoeffun(location,state)
    N = 4;
    nr = length(location.x);
    f = zeros(N,nr);

    D_f = 1.96*(10^(-6)); % Lit
    D_m = 9*(10^(-6)); % Lit
    D_p = 9*(10^(-6)); % Lit
    delta_c = 0.00001; % Estimate
    delta_f = 0.000139; % Lit
    delta_m = log(2)/49; % Lit
    delta_p = log(2)/4; % Lit
    
    sigma_c = 0.0183; % Lit
    c_0 = 0.68; % Lit
    f_max = 15; % Estimate
    sigma_f = 0.15;% Estimate
    f_0 = 0.632; % Lit
    k_p = 0.235; % % Lit

    eta = 9.144*10^(6); % Lit
    mu = 3.6*(10)^7; % Lit
    k_m = 2.523*10^(-10); % Lit
    m_thresh = 9*10^(-10); % Estimate
    
    sigma_C = f_0/c_0;
    eta_bar = eta*m_thresh/sigma_c;
    delta_C = delta_c/sigma_c;
    sigma_F = sigma_f/sigma_c;
    delta_F = delta_f/sigma_c;
    k_M = k_m/(m_thresh*sigma_c);
    k_M = 110;
    mu_bar = mu*m_thresh/sigma_c;
    delta_M = delta_m/sigma_c;
    k_P = k_p/sigma_c;
    delta_P = delta_p/sigma_c;

    f(1,:) = sigma_C.*state.u(2,:).*(1-state.u(1,:)).*(1+(f_max.*state.u(3,:).*exp(1-state.u(3,:))))...
        - eta_bar.*state.u(1,:).*state.u(3,:) - delta_C.*state.u(1,:);

    f(2,:) = sigma_F.*state.u(2,:).*(1-state.u(2,:)) - delta_F.*state.u(2,:);

    f(3,:) = k_M.*state.u(2,:).*(1-state.u(1,:))...
        - mu_bar.*state.u(3,:).*state.u(4,:) - delta_M.*state.u(3,:);

    f(4,:) = k_P.*state.u(2,:).*state.u(3,:)...
        - mu_bar.*state.u(3,:).*state.u(4,:) - delta_P.*state.u(4,:);
end

 



