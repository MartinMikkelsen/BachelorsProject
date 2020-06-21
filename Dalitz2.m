%Dalitz 2
clc; clear all; close all
M  = 3749.41;                                   %MeV
m1 = 3728.4;                                   %MeV
m2 = 0.5109989461;                              %MeV
m3 = 0.5109989461;                              %MeV
lambda = @(x,y,z) (x-y-z).^2-4*y*z;             %Triangle function
s = M^2;                                        %Mass of the decaying particle squared
s1 = linspace((m1+m2).^2,(sqrt(s)-m3)^2,1000); %s1 = s12
s2 = linspace((m2+m3).^2,(sqrt(s)-m1)^2,1000); %s2 = s23
s3 = linspace((m1+m3).^2,(sqrt(s)-m2)^2,1000); %s3 = s31
C = linspecer(4); 

s11 = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)-lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)); %upper half of boundary curve
s12 = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)+lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)); %lower half of boundary curve
hold on
plot(s2,s11(s2),'r','LineWidth',5,'color',C(1,:))
plot(s2,s12(s2),'r','LineWidth',5,'color',C(1,:))

xlabel('$s_{2}[MeV^2]$ ','Interpreter','latex', 'Fontsize',42);
ylabel('$s_{1} [MeV^2]$ ','Interpreter','latex','Fontsize',42);
%legend("Boundary curve for $D_s^+ \rightarrow \pi^+ +\pi^- + \pi^+$",'Interpreter','latex','Fontsize',42,'Location', 'Best')

set(gca,'FontSize',42)
%title('Dalitz plot','Interpreter','latex', 'Fontsize',40)
ax.LineWidth = 2;
ax.FontSize = 42;
ax.TickLabelInterpreter = 'LaTeX';
set(gcf, 'PaperPositionMode', 'auto');
%print -depsc2 Dalitzplot2.eps
%% Now calculate the projection onto s1,s2,s3 for spin = 0,1,3. This is split into 3 parts for 3 different angular dist.
mr = 1/3*0.0167*1000;                                %Pole mass in MeV
mr1 = 1/3*mr;
mr2 = 1.2*mr;
gamma01 =   2;                         %Decay width in MeV
A1 =@(a) sqrt(a)./(mr^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner
%A4 =@(a) sqrt(a)./(mr1^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner
%A5 =@(a) sqrt(a)./(mr2^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner

eta2 = ((mr^2+s1-m1)./(2*mr^2)).^2-1;           %Relativistic correction
%Decay through channel 1
%% spin = 0 
T0 = 1;
Matrixelement = @(b) T0.*A1(b);
Newmatrixelement1 = @(b) T0.*A4(b);
Newmatrixelement2 = @(b) T0.*A5(b);
hold on

plot(s2,Matrixelement(s2).*conj(Matrixelement(s2)),'color',C(1,:))
%plot(s2,Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)),'color',C(2,:))
%plot(s2,Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)),'color',C(3,:))
%%
costheta = @(s1) 2.*(sqrt(s1)-min(sqrt(s1)))./(max(sqrt(s1))-min(sqrt(s1)))-1;
lpol = @(s1) legendreP(1,costheta(s1));
hold on
theta = linspace(0,pi,1000);
plot(theta,lpol(s1).^2)
A1 =@(a) sqrt(a)./(mr^2-a-1i*gamma01*sqrt(a));
camroll(-90)
%%
[x,y] = meshgrid(s2,[s11(s2),s12(s2)]);
hold on


mscale1 = @(x) (Matrixelement(s2) - min(Matrixelement(s2)))./( max(Matrixelement(s2)) - min(Matrixelement(s2))); %normalize
z = (mscale1(x).*conj(mscale1(x))).*abs(costheta(y));%*abs(1/2.*(3.*costheta(y).^2-1)); %remember spin

bad=(y>=s11(s2)) | (y<=s12(s2));
z(bad)=nan;

znew = z-min(z)/(max(z)-min(z));

contourf(x,real(y),real(z),'edgecolor','none')
plot(s2,s11(s2),'r','LineWidth',2)
plot(s2,s12(s2),'r','LineWidth',2)
colormap(linspecer);
%%
a1 = subplot(3,4,[1 2 3]);
hold on
box on
plot(s2,(Matrixelement(s2).*conj(Matrixelement(s2))-min(Matrixelement(s2).*conj(Matrixelement(s2))))./(max(Matrixelement(s2).*conj(Matrixelement(s2)))-min(Matrixelement(s2).*conj(Matrixelement(s2)))),'LineWidth',3,'color',C(1,:))
%plot(s2,(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2))-min(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2))))./(max(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)))-min(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)))),'--','LineWidth',3,'color',C(2,:))
%plot(s2,(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2))-min(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2))))./(max(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)))-min(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)))),'--','LineWidth',3,'color',C(3,:))
%legend("$m_X$","$1/3 m_X$","$1.2 m_X$",'Interpreter','latex','Fontsize',24,'Location', 'Best')
a1.TickLabelInterpreter = 'LaTeX';
set(gca,'xtick',[])
set(gca,'FontSize',50)
a1.LineWidth = 4;
ylabel('Events','Interpreter','latex');
a2 = subplot(3,4,[8 12]);
plot(theta,lpol(s1).^2,'LineWidth',3,'color',C(1,:))
box on
a2.TickLabelInterpreter = 'LaTeX';
a2.LineWidth = 4;
camroll(-90)
set(gca,'xtick',[])
%set(gca,'ytick',[0.2 0.4 0.6 0.8 1])
set(gca,'FontSize',50)
ylabel('Events','Interpreter','latex', 'Fontsize',50);
a3 = subplot(3,4,[5 6 7 9 10 11]);
a3.FontSize = 50;
hold on
contourf(x,real(y),real(z),'edgecolor','none')
plot(s2,s11(s2),'k','LineWidth',2)
plot(s2,s12(s2),'k','LineWidth',2)
a3.FontSize = 42;
a2.LineWidth = 4;
a3.TickLabelInterpreter = 'LaTeX';
set(gca,'FontSize',50)
set(gcf, 'PaperPositionMode', 'auto');
xlabel('$s_{2} [MeV^2]$','Interpreter','latex', 'Fontsize',50,'LineWidth',9);
ylabel('$s_{1} [MeV^2]$','Interpreter','latex','Fontsize',50);
a3.LineWidth = 4;
bad=(y>=s11(s2)) | (y<=s12(s2));
z(bad)=nan;
colormap(linspecer);
set(gca,'TickLabelInterpreter','latex')
box on
h = colorbar;
h.TickLabelInterpreter = 'LaTeX';


%legend("Boundary curve for $He^* \rightarrow He+e^+ + e^-$",'Interpreter','latex','Fontsize',15)
%%
E2 = (s+m2^2-s3)./(2*sqrt(s));
E3 = (s+m3^2-s1)./(2*sqrt(s));
P2 = (lambda(s,m2^2,s3).^(1/2))/2.*sqrt(s);
P3 = (lambda(s,m3^2,s1).^(1/2))/2.*sqrt(s);
s2theta140 = m2^2+m3^2+2.*E2.*E3-2.*P2.*P3.*cosd(140)
%%
costheta = @(s1) 2.*(sqrt(s1)-min(sqrt(s1)))./(max(sqrt(s1))-min(sqrt(s1)))-1;
lpol = @(s1) legendreP(0,costheta(s1));
hold on
box on
ax = gca;
lpol1 = @(s1) legendreP(1,costheta(s1));
lpol2 = @(s1) legendreP(2,costheta(s1));
lpol3 = @(s1) legendreP(3,costheta(s1));
theta = linspace(pi,0,1000);
b1 = plot(theta,lpol(s1).^2,'r','LineWidth',6,'color',C(1,:));
b2 = plot(theta,lpol1(s1).^2,'--m','LineWidth',6,'color',C(2,:));
b3 = plot(theta,lpol2(s1).^2,'b','LineWidth',6,'color',C(3,:));
b4 = plot(theta,lpol3(s1).^2,'--k','LineWidth',6,'color',C(4,:));
ylim([0 1.1])
xlim([0 pi])
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'-1','0','1'})
A1 =@(a) sqrt(a)./(mr^2-a-1i*gamma01*sqrt(a));
ax.LineWidth = 3;
ax.FontSize = 42;
set(gca,'FontSize',32)
ax.TickLabelInterpreter = 'LaTeX';
set(gcf, 'PaperPositionMode', 'auto');
xlabel('$\cos \theta$','Fontsize',123,'Interpreter','latex');
ylabel('Arb. units','Interpreter','latex','Fontsize',32);
legend("Spin 0","Spin 1","Spin 2","Spin 3",'Interpreter','latex','Fontsize',42,'Location', 'Best')
set(gca,'TickLabelInterpreter','latex')
ax.FontSize = 42;
ax.TickLabelInterpreter = 'LaTeX';
set(gca,'FontSize',60,'LineWidth',5)
set(gcf, 'PaperPositionMode', 'auto');

%%
[xrect,yrect] = meshgrid(s2,[s11(s2),s12(s2)]);
 
[x,y]=deal(xrect,yrect);
bad=(y>=s11(s2)) | (y<=s12(s2));

x(bad)=[];
y(bad)=[];
mscale1 = @(s2) (Matrixelement(s2) - min(Matrixelement(s2)))./( max(Matrixelement(s2)) - min(Matrixelement(s2))); %normalize
z = (mscale1(x).^2).*abs(1/2.*(3.*costheta(y).^2-1)); %remember spin

zrect=nan(size(xrect));
zrect(~bad)=z;

contourf(xrect,real(yrect),real(zrect),'edgecolor','none')
colormap(linspecer);

box on
colorbar