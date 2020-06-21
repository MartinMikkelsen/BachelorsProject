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
C = linspecer(5); 

s11 = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)-lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)); %upper half of boundary curve
s12 = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)+lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)); %lower half of boundary curve

s11test = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)-cosd(25).*((lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)))); %upper half of boundary curve
s12test = @(s2) m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)+cosd(25).*((lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2)))); %lower half of boundary curve
hold on
plot(s2,s11(s2),'r','LineWidth',5,'color',C(1,:))
plot(s2,s12(s2),'r','LineWidth',5,'color',C(1,:))

plot(s2,s11test(s2),'r','LineWidth',5,'color',C(1,:))
plot(s2,s12test(s2),'r','LineWidth',5,'color',C(1,:))



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
mr = 0.0167*1000;                                %Pole mass in MeV
mr1 = 1/3*mr;
mr2 = 1.25*mr;
gamma01 =   2;                         %Decay width in MeV
A1 =@(a) sqrt(a)./(mr^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner
A4 =@(a) sqrt(a)./(mr1^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner
A5 =@(a) sqrt(a)./(mr2^2-a-1i*gamma01*sqrt(a));  %Breit-Wigner

eta2 = ((mr^2+s1-m1)./(2*mr^2)).^2-1;           %Relativistic correction
%Decay through channel 1
%% spin = 0 
T0 = 1;
Matrixelement = @(b) T0.*A1(b);
Newmatrixelement1 = @(b) T0.*A4(b);
Newmatrixelement2 = @(b) T0.*A5(b);
hold on

plot(s2,Matrixelement(s2).*conj(Matrixelement(s2)),'color',C(1,:))
plot(s2,Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)),'color',C(2,:))
plot(s2,Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)),'color',C(3,:))
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
z = (mscale1(x).*conj(mscale1(x))).*abs(costheta(y));%abs(1/2.*(3.*costheta(y).^2-1)); %remember spin

bad=(y>=s11(s2)) | (y<=s12(s2));
z(bad)=nan;

znew = z-min(z)/(max(z)-min(z));

contourf(x,real(y),real(z),'edgecolor','none')
plot(s2,s11(s2),'r','LineWidth',2)
plot(s2,s12(s2),'r','LineWidth',2)
colormap(linspecer);
%%
a1 = subplot(3,4,[1 2 3])
hold on
box on
a1.LineWidth = 3;
a1.TickLabelInterpreter = 'LaTeX';
plot(s2,(Matrixelement(s2).*conj(Matrixelement(s2))-min(Matrixelement(s2).*conj(Matrixelement(s2))))./(max(Matrixelement(s2).*conj(Matrixelement(s2)))-min(Matrixelement(s2).*conj(Matrixelement(s2)))),'LineWidth',3,'color',C(1,:))
plot(s2,(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2))-min(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2))))./(max(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)))-min(Newmatrixelement1(s2).*conj(Newmatrixelement1(s2)))),'--','LineWidth',3,'color',C(2,:))
plot(s2,(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2))-min(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2))))./(max(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)))-min(Newmatrixelement2(s2).*conj(Newmatrixelement2(s2)))),'--','LineWidth',3,'color',C(3,:))
legend("$m_X$","$1/3 m_X$","$1.2 m_X$",'Interpreter','latex','Fontsize',24,'Location', 'Best')
lgd.FontSize = 24;
%ax.FontSize = 40;
set(gca,'xtick',[])
set(gca,'FontSize',24)
%set(gca,'ytick',[])
ylabel('Events','Interpreter','latex', 'Fontsize',40);
a2 =subplot(3,4,[8 12])
plot(theta,lpol(s1).^2,'LineWidth',3,'color',C(1,:))
box on
a2.LineWidth = 3;
a2.TickLabelInterpreter = 'LaTeX';
camroll(-90)
set(gca,'xtick',[])
%set(gca,'ytick',[0.2 0.4 0.6 0.8 1])
set(gca,'FontSize',24)
ylabel('Events','Interpreter','latex', 'Fontsize',40);
set(gca,'TickLabelInterpreter','latex')
a3 =subplot(3,4,[5 6 7 9 10 11])
%ax.FontSize = 40;
hold on
contourf(x,real(y),real(z),'edgecolor','none')
plot(s2,s11(s2),'k','LineWidth',2)
plot(s2,s12(s2),'k','LineWidth',2)
%q1 = plot(s2,s11test(s2),'--','LineWidth',5,'color',C(4,:))
%q2 = plot(s2,s12test(s2),'--','LineWidth',5,'color',C(4,:))
a3.LineWidth = 3;
%ax.FontSize = 30;

%ax.TickLabelInterpreter = 'LaTeX';
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',30)
set(gcf, 'PaperPositionMode', 'auto');
xlabel('$s_{2} [MeV^2]$','Interpreter','latex', 'Fontsize',40);
ylabel('$s_{1} [MeV^2]$','Interpreter','latex','Fontsize',40);
bad=(y>=s11(s2)) | (y<=s12(s2));
z(bad)=nan;
colormap(linspecer);

box on
colorbar
%legend([q1 q2],"$\theta^{CM}=25^\circ$",'Interpreter','latex','Fontsize',24,'Location', 'Best')

%legend("Boundary curve for $He^* \rightarrow He+e^+ + e^-$",'Interpreter','latex','Fontsize',15)
%ax = gca;
%exportgraphics(ax,'HeJ2.pdf','ContentType','vector')

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

%%

% costheta12 = ((s+m1^2-s2)*(s+m2^2-s3)+2*s*(m1^2+m2^2-s1))/sqrt((lambda(s,m1^2,s2))*sqrt(lambda(s,m2^2,s3)))
% angle(costheta12)
% 
% costheta23 = ((s+m2^2-s3)*(s+m3^2-s1)+2*s*(m2^2+m3^2-s2))/sqrt((lambda(s,m1^2,s3))*sqrt(lambda(s,m3^2,s1)))
% angle(costheta23)
% 
% costheta31 = ((s+m3^2-s1)*(s+m1^2-s2)+2*s*(m1^2+m3^2-s3))/sqrt((lambda(s,m1^2,s1))*sqrt(lambda(s,m1^2,s2)))
% angle(costheta31)
%%
format short
s1 = 1.398*10^7;
s2 = 220;
s3 = s+m1^2+m2^2+m3^2-s1-s2;
E1 = (s+m1^2-s2)/(2*sqrt(s));
E2 = (s+m2^2-s3)/(2*sqrt(s));
E3 = (s+m3^2-s1)/(2*sqrt(s));
P1 = (sqrt(lambda(s,m1^2,s2)))/(2*sqrt(s));
P2 = (sqrt(lambda(s,m2^2,s3)))/(2*sqrt(s));
P3 = (sqrt(lambda(s,m3^2,s1)))/(2*sqrt(s));

costheta12test = (m1^2+m2^2-s1+2*E1*E2)/(2*P1*P2);
costheta23test = (m2^2+m3^2-s2+2*E2*E3)/(2*P2*P3);
costheta31test = (m1^2+m3^2-s3+2*E1*E3)/(2*P1*P3);

thetas1 = acosd(costheta12test)
thetas2 = acosd(costheta23test)
thetas3 = acosd(costheta31test)

test = m1^2+m2^2-1./(2*s2).*((s2-s+m1^2).*(s2+m2^2-m3^2)-cosd(105).*((lambda(s2,s,m1^2).^(1/2).*lambda(s2,m2^2,m3^2).^(1/2))))
%%
A = [180  1828.5310346476
168.52165904547  1828.9420677382
163.73979529169  1829.3531008289
160.0515564112  1829.7641339196
156.92608193437  1830.1751670103
154.15806723683  1830.586200101
151.64236342367  1830.9972331917
149.31658289102  1831.4082662824
147.14011962111  1831.8192993731
145.08479375256  1832.2303324638
143.13010235416  1832.6413655545
141.26057540214  1833.0523986452
139.46419788868  1833.4634317359
137.73141557043  1833.8744648265
136.05448043769  1834.2854979172
134.42700400081  1834.6965310079
132.8436430436  1835.1075640986
131.29987279171  1835.5185971893
129.79181949956  1835.92963028
128.31613447367  1836.3406633707
126.86989764584  1836.7516964614
125.45054263918  1837.1627295521
124.05579774257  1837.5737626428
122.68363884626  1837.9847957335
121.33225149759  1838.3958288242
120  1838.8068619149
118.68540201412  1839.2178950055
117.38710750265  1839.6289280962
116.10388113734  1840.0399611869
114.8345874897  1840.4509942776
113.5781784782  1840.8620273683
112.33368265781  1841.273060459
111.10019602409  1841.6840935497
109.87687407008  1842.0951266404
108.66292488494  1842.5061597311
107.45760312372  1842.9171928218
106.26020470831  1843.3282259125
105.07006214489  1843.7392590032
103.88654036263  1844.1502920938
102.7090329944  1844.5613251845
101.53695903282  1844.9723582752
100.36975980548  1845.3833913659
99.206896221346  1845.7944244566
98.047846247312  1846.2054575473
96.892102579346  1846.616490638
95.739170477267  1847.0275237287
94.588565735786  1847.4385568194
93.439812767515  1847.8495899101
92.292442775956  1848.2606230008
91.145991998389  1848.6716560915
90  1849.0826891821
88.854008001611  1849.4937222728
87.707557224044  1849.9047553635
86.560187232485  1850.3157884542
85.411434264214  1850.7268215449
84.260829522733  1851.1378546356
83.107897420654  1851.5488877263
81.952153752688  1851.959920817
80.793103778654  1852.3709539077
79.630240194523  1852.7819869984
78.463040967185  1853.1930200891
77.290967005605  1853.6040531798
76.113459637371  1854.0150862705
74.929937855111  1854.4261193611
73.739795291688  1854.8371524518
72.542396876278  1855.2481855425
71.337075115058  1855.6592186332
70.123125929921  1856.0702517239
68.899803975907  1856.4812848146
67.666317342195  1856.8923179053
66.421821521798  1857.303350996
65.165412510298  1857.7143840867
63.89611886266  1858.1254171774
62.612892497346  1858.5364502681
61.314597985881  1858.9474833588
60  1859.3585164494
58.667748502406  1859.7695495401
57.316361153742  1860.1805826308
55.944202257432  1860.5916157215
54.549457360825  1861.0026488122
53.130102354156  1861.4136819029
51.683865526334  1861.8247149936
50.208180500443  1862.2357480843
48.700127208294  1862.646781175
47.156356956404  1863.0578142657
45.572995999194  1863.4688473564
43.945519562309  1863.8798804471
42.268584429572  1864.2909135378
40.535802111317  1864.7019466284
38.739424597856  1865.1129797191
36.869897645844  1865.5240128098
34.915206247444  1865.9350459005
32.859880378889  1866.3460789912
30.683417108976  1866.7571120819
28.357636576328  1867.1681451726
25.841932763167  1867.5791782633
23.073918065631  1867.990211354
19.948443588803  1868.4012444447
16.260204708312  1868.8122775354
11.478340954534  1869.2233106261
0  1869.6343437167];
A1 = A(:,1);
A2 = A(:,2);
plot(A1,A2,'r','LineWidth',5,'color',C(1,:))
xlabel('$\theta_{cm} $[deg] ','Interpreter','latex', 'Fontsize',42);
ylabel('$E_3 $[MeV]','Interpreter','latex','Fontsize',42);
legend("$p+^3H \rightarrow ^4He^*+X, E_X(4)=21.01 MeV, E_k(p)=400 KeV$",'Interpreter','latex','Fontsize',42,'Location', 'Best')
set(gca,'FontSize',42)
ax.LineWidth = 2;
ax.FontSize = 42;
ax.TickLabelInterpreter = 'LaTeX';
set(gcf, 'PaperPositionMode', 'auto');
