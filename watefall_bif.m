%% Bifurcation diagrams obtained using XPP at the below specified values
data1 = load('allinfo_kP_0_1.dat');
data2 = load('allinfo_kP_0_2.dat');
data3 = load('allinfo_kP_0_3.dat');
data4 = load('allinfo_kP_0_4.dat');
data5 = load('allinfo_kP_0_5.dat');
data6 = load('allinfo_kP_0_6.dat');
data7 = load('allinfo_kP_0_7.dat');
data8 = load('allinfo_kP_0_8.dat');
data9 = load('allinfo_kP_0_9.dat');
data10 = load('allinfo_kP_1_0.dat');
data11 = load('allinfo_kP_1_1.dat');
data12 = load('allinfo_kP_1_2.dat');
data13 = load('allinfo_kP_1_3.dat');
data14 = load('allinfo_kP_1_4.dat');
data15 = load('allinfo_kP_1_5.dat');
data16 = load('allinfo_kP_1_6.dat');
data17 = load('allinfo_kP_1_7.dat');
data18 = load('allinfo_kP_1_8.dat');
data19 = load('allinfo_kP_1_9.dat');
data20 = load('allinfo_kP_2_0.dat');
data21 = load('allinfo_kP_2_1.dat');
data22 = load('allinfo_kP_2_2.dat');
data23 = load('allinfo_kP_2_3.dat');
data24 = load('allinfo_kP_2_4.dat');
data25 = load('allinfo_kP_2_5.dat');

xval1 = data1(:,4);
xval2 = data2(:,4);
xval3 = data3(:,4);
xval4 = data4(:,4);
xval5 = data5(:,4);
xval6 = data6(:,4);
xval7 = data7(:,4);
xval8 = data8(:,4);
xval9 = data9(:,4);
xval10 = data10(:,4);
xval11 = data11(:,4);
xval12 = data12(:,4);
xval13 = data13(:,4);
xval14 = data14(:,4);
xval15 = data15(:,4);
xval16 = data16(:,4);
xval17 = data17(:,4);
xval18 = data18(:,4);
xval19 = data19(:,4);
xval20 = data20(:,4);
xval21 = data21(:,4);
xval22 = data22(:,4);
xval23 = data23(:,4);
xval24 = data24(:,4);
xval25 = data25(:,4);

yval1 = data1(:,7);
yval2 = data2(:,7);
yval3 = data3(:,7);
yval4 = data4(:,7);
yval5 = data5(:,7);
yval6 = data6(:,7);
yval7 = data7(:,7);
yval8 = data8(:,7);
yval9 = data9(:,7);
yval10 = data10(:,7);
yval11 = data11(:,7);
yval12 = data12(:,7);
yval13 = data13(:,7);
yval14 = data14(:,7);
yval15 = data15(:,7);
yval16 = data16(:,7);
yval17 = data17(:,7);
yval18 = data18(:,7);
yval19 = data19(:,7);
yval20 = data20(:,7);
yval21 = data21(:,7);
yval22 = data22(:,7);
yval23 = data23(:,7);
yval24 = data24(:,7);
yval25 = data25(:,7);

M1 = [xval1,yval1];
M2 = [xval2,yval2];
M3 = [xval3,yval3];
M4 = [xval4,yval4];
M5 = [xval5,yval5];
M6 = [xval6,yval6];
M7 = [xval7,yval7];
M8 = [xval8,yval8];
M9 = [xval9,yval9];
M10 = [xval10,yval10];
M11 = [xval11,yval11];
M12 = [xval12,yval12];
M13 = [xval13,yval13];
M14 = [xval14,yval14];
M15 = [xval15,yval15];
M16 = [xval16,yval16];
M17 = [xval17,yval17];
M18 = [xval18,yval18];
M19 = [xval19,yval19];
M20 = [xval20,yval20];
M21 = [xval21,yval21];
M22 = [xval22,yval22];
M23 = [xval23,yval23];
M24 = [xval24,yval24];
M25 = [xval25,yval25];

% X = [M1(:,1),M2(:,1),M3(:,1),M4(:,1),M5(:,1),M6(:,1),M7(:,1),M8(:,1),...
%     M9(:,1),M10(:,1),M11(:,1),M12(:,1),M13(:,1),M14(:,1),M15(:,1),M16(:,1),...
%     M17(:,1),M18(:,1),M19(:,1),M20(:,1),M21(:,1),M22(:,1),M23(:,1),M24(:,1),...
%     M24(:,1)];
% X = transpose(X);
% 
% Y = [M1(:,2),M2(:,2),M3(:,2),M4(:,2),M5(:,2),M6(:,2),M7(:,2),M8(:,2),...
%     M9(:,2),M10(:,2),M11(:,2),M12(:,2),M13(:,2),M14(:,2),M15(:,2),M16(:,2),...
%     M17(:,2),M18(:,2),M19(:,2),M20(:,2),M21(:,2),M22(:,2),M23(:,2),M24(:,2),...
%     M24(:,2)];
% Y = transpose(Y);
% 
% Z = [0.1*ones(200,1),0.2*ones(200,1),0.3*ones(200,1),0.4*ones(200,1),...
%     0.5*ones(200,1),0.6*ones(200,1),0.7*ones(200,1),0.8*ones(200,1),...
%     0.9*ones(200,1),1.0*ones(200,1),1.1*ones(200,1),1.2*ones(200,1),...
%     1.3*ones(200,1),1.4*ones(200,1),1.5*ones(200,1),1.6*ones(200,1),...
%     1.7*ones(200,1),1.8*ones(200,1),1.9*ones(200,1),2.0*ones(200,1),...
%     2.1*ones(200,1),2.2*ones(200,1),2.3*ones(200,1),2.4*ones(200,1),...
%     2.5*ones(200,1)];
% Z = transpose(Z);

% figure(2)
% waterfall(X,Z,Y)

U1 = nan(length(data1),width(data1));
for i = 1:length(data1)
    if data1(i,1) == 2
        U1(i,:) = data1(i,:);
    end
end
hold on
plot3(data1(:,4),0.1*ones(200,1),data1(:,7),'color','k','linewidth',1.5)
plot3(U1(:,4),0.1*ones(200,1),U1(:,7),'color','r','linewidth',1.5)

U2 = nan(length(data2),width(data2));
for i = 1:length(data2)
    if data2(i,1) == 2
        U2(i,:) = data2(i,:);
    end
end
hold on
plot3(data2(:,4),0.2*ones(200,1),data2(:,7),'color','k','linewidth',1.5)
plot3(U2(:,4),0.2*ones(200,1),U2(:,7),'color','r','linewidth',1.5)

U3 = nan(length(data3),width(data3));
for i = 1:length(data3)
    if data3(i,1) == 2
        U3(i,:) = data3(i,:);
    end
end
hold on
plot3(data3(:,4),0.3*ones(200,1),data3(:,7),'color','k','linewidth',1.5)
plot3(U3(:,4),0.3*ones(200,1),U3(:,7),'color','r','linewidth',1.5)

U4 = nan(length(data4),width(data4));
for i = 1:length(data4)
    if data4(i,1) == 2
        U4(i,:) = data4(i,:);
    end
end
hold on
plot3(data4(:,4),0.4*ones(200,1),data4(:,7),'color','k','linewidth',1.5)
plot3(U4(:,4),0.4*ones(200,1),U4(:,7),'color','r','linewidth',1.5)

U5 = nan(length(data5),width(data5));
for i = 1:length(data5)
    if data5(i,1) == 2
        U5(i,:) = data5(i,:);
    end
end
hold on
plot3(data5(:,4),0.5*ones(200,1),data5(:,7),'color','k','linewidth',1.5)
plot3(U5(:,4),0.5*ones(200,1),U5(:,7),'color','r','linewidth',1.5)

U6 = nan(length(data6),width(data6));
for i = 1:length(data6)
    if data6(i,1) == 2
        U6(i,:) = data6(i,:);
    end
end
hold on
plot3(data6(:,4),0.6*ones(200,1),data6(:,7),'color','k','linewidth',1.5)
plot3(U6(:,4),0.6*ones(200,1),U6(:,7),'color','r','linewidth',1.5)


U7 = nan(length(data7),width(data7));
for i = 1:length(data7)
    if data7(i,1) == 2
        U7(i,:) = data7(i,:);
    end
end
hold on
plot3(data7(:,4),0.7*ones(200,1),data7(:,7),'color','k','linewidth',1.5)
plot3(U7(:,4),0.7*ones(200,1),U7(:,7),'color','r','linewidth',1.5)

U8 = nan(length(data8),width(data8));
for i = 1:length(data8)
    if data8(i,1) == 2
        U8(i,:) = data8(i,:);
    end
end
hold on
plot3(data8(:,4),0.8*ones(200,1),data8(:,7),'color','k','linewidth',1.5)
plot3(U8(:,4),0.8*ones(200,1),U8(:,7),'color','r','linewidth',1.5)

U9 = nan(length(data9),width(data9));
for i = 1:length(data9)
    if data9(i,1) == 2
        U9(i,:) = data9(i,:);
    end
end
hold on
plot3(data9(:,4),0.9*ones(200,1),data9(:,7),'color','k','linewidth',1.5)
plot3(U9(:,4),0.9*ones(200,1),U9(:,7),'color','r','linewidth',1.5)

U10 = nan(length(data10),width(data10));
for i = 1:length(data10)
    if data10(i,1) == 2
        U10(i,:) = data10(i,:);
    end
end
hold on
plot3(data10(:,4),1*ones(200,1),data10(:,7),'color','k','linewidth',1.5)
plot3(U10(:,4),1*ones(200,1),U10(:,7),'color','r','linewidth',1.5)

U11 = nan(length(data11),width(data11));
for i = 1:length(data11)
    if data11(i,1) == 2
        U11(i,:) = data11(i,:);
    end
end
hold on
plot3(data11(:,4),1.1*ones(200,1),data11(:,7),'color','k','linewidth',1.5)
plot3(U11(:,4),1.1*ones(200,1),U11(:,7),'color','r','linewidth',1.5)

U12 = nan(length(data12),width(data12));
for i = 1:length(data12)
    if data12(i,1) == 2
        U12(i,:) = data12(i,:);
    end
end
hold on
plot3(data12(:,4),1.2*ones(200,1),data12(:,7),'color','k','linewidth',1.5)
plot3(U12(:,4),1.2*ones(200,1),U12(:,7),'color','r','linewidth',1.5)

U13 = nan(length(data13),width(data13));
for i = 1:length(data13)
    if data13(i,1) == 2
        U13(i,:) = data13(i,:);
    end
end
hold on
plot3(data13(:,4),1.3*ones(200,1),data13(:,7),'color','k','linewidth',1.5)
plot3(U13(:,4),1.3*ones(200,1),U13(:,7),'color','r','linewidth',1.5)

U14 = nan(length(data14),width(data14));
for i = 1:length(data14)
    if data14(i,1) == 2
        U14(i,:) = data14(i,:);
    end
end
hold on
plot3(data14(:,4),1.4*ones(200,1),data14(:,7),'color','k','linewidth',1.5)
plot3(U14(:,4),1.4*ones(200,1),U14(:,7),'color','r','linewidth',1.5)

U15 = nan(length(data15),width(data15));
for i = 1:length(data15)
    if data15(i,1) == 2
        U15(i,:) = data15(i,:);
    end
end
hold on
plot3(data15(:,4),1.5*ones(200,1),data15(:,7),'color','k','linewidth',1.5)
plot3(U15(:,4),1.5*ones(200,1),U15(:,7),'color','r','linewidth',1.5)

U16 = nan(length(data16),width(data16));
for i = 1:length(data16)
    if data16(i,1) == 2
        U16(i,:) = data16(i,:);
    end
end
hold on
plot3(data16(:,4),1.6*ones(200,1),data16(:,7),'color','k','linewidth',1.5)
plot3(U16(:,4),1.6*ones(200,1),U16(:,7),'color','r','linewidth',1.5)


U17 = nan(length(data17),width(data17));
for i = 1:length(data17)
    if data17(i,1) == 2
        U17(i,:) = data17(i,:);
    end
end
hold on
plot3(data17(:,4),1.7*ones(135,1),data17(:,7),'color','k','linewidth',1.5)
plot3(U17(:,4),1.7*ones(135,1),U17(:,7),'color','r','linewidth',1.5)

U18 = nan(length(data18),width(data18));
for i = 1:length(data18)
    if data18(i,1) == 2
        U18(i,:) = data18(i,:);
    end
end
hold on
plot3(data18(:,4),1.8*ones(200,1),data18(:,7),'color','k','linewidth',1.5)
plot3(U18(:,4),1.8*ones(200,1),U18(:,7),'color','r','linewidth',1.5)

U19 = nan(length(data19),width(data19));
for i = 1:length(data19)
    if data19(i,1) == 2
        U19(i,:) = data19(i,:);
    end
end
hold on
plot3(data19(:,4),1.9*ones(200,1),data19(:,7),'color','k','linewidth',1.5)
plot3(U19(:,4),1.9*ones(200,1),U19(:,7),'color','r','linewidth',1.5)

U20 = nan(length(data20),width(data20));
for i = 1:length(data20)
    if data20(i,1) == 2
        U20(i,:) = data20(i,:);
    end
end
hold on
plot3(data20(:,4),2*ones(200,1),data20(:,7),'color','k','linewidth',1.5)
plot3(U20(:,4),2*ones(200,1),U20(:,7),'color','r','linewidth',1.5)

U21 = nan(length(data21),width(data21));
for i = 1:length(data21)
    if data21(i,1) == 2
        U21(i,:) = data21(i,:);
    end
end
hold on
plot3(data21(:,4),2.1*ones(200,1),data21(:,7),'color','k','linewidth',1.5)
plot3(U21(:,4),2.1*ones(200,1),U21(:,7),'color','r','linewidth',1.5)

U22 = nan(length(data22),width(data22));
for i = 1:length(data22)
    if data22(i,1) == 2
        U22(i,:) = data22(i,:);
    end
end
hold on
plot3(data22(:,4),2.2*ones(200,1),data22(:,7),'color','k','linewidth',1.5)
plot3(U22(:,4),2.2*ones(200,1),U22(:,7),'color','r','linewidth',1.5)

U23 = nan(length(data23),width(data23));
for i = 1:length(data23)
    if data23(i,1) == 2
        U23(i,:) = data23(i,:);
    end
end
hold on
plot3(data23(:,4),2.3*ones(200,1),data23(:,7),'color','k','linewidth',1.5)
plot3(U23(:,4),2.3*ones(200,1),U23(:,7),'color','r','linewidth',1.5)

U24 = nan(length(data24),width(data24));
for i = 1:length(data24)
    if data24(i,1) == 2
        U24(i,:) = data24(i,:);
    end
end
hold on
plot3(data24(:,4),2.4*ones(200,1),data24(:,7),'color','k','linewidth',1.5)
plot3(U24(:,4),2.4*ones(200,1),U24(:,7),'color','r','linewidth',1.5)

U25 = nan(length(data25),width(data25));
for i = 1:length(data25)
    if data25(i,1) == 2
        U25(i,:) = data25(i,:);
    end
end
hold on
plot3(data25(:,4),2.5*ones(200,1),data25(:,7),'color','k','linewidth',1.5)
plot3(U25(:,4),2.5*ones(200,1),U25(:,7),'color','r','linewidth',1.5)

grid on 
set(gca,'fontsize',16)
xlabel('$\delta_C$','FontSize',18,'FontWeight','bold','interpreter','latex')
ylabel('$k_P$','FontSize',18,'FontWeight','bold','interpreter','latex')
zlabel('$\overline{C}$','FontSize',18,'FontWeight','bold','interpreter','latex')
xlim([0,50])
