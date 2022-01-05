% Exemplary call of QPstab to investigate the stability of the origin of
% the Matthieu equation

% q_dd+c*q_d+[k+eps*(ampls(1)*cos(Om_1 *t)+ampls(1)*cos(Om_2 *t) ) ] q=0


% First order form of the equation of variation must be as separate function
% in a m-file (e.g. eq_of_var.m for the above example)
clear all
close all



irr=sqrt(3)*10^-6;
ws=linspace(irr,1.5+irr,100);
ks=linspace(0,0.5,100);

% select parameters
eps=0.1;
ampls=[1   1];
c=0.01;

% Maximal number of periodis for which we integrate the equations of 
% variation, i.e. Tmax=Nmax*T=Nmax*2*pi/Oms(1)
Nmax=400;
% estimated number of neighbourhoods along each edge of the cubic section
% of the trap, we set N_edge=10 in all our calculations
N_edge=10;
% set tolerance, if estimated neighbourhood, i.e. Delta(phi_0) in eq. (20),
% shrinks below delta_tol the calculation are terminated 
delta_tol= 10^-5;


% Dimensions of the phase space of eq. (1), for the first order equivalent
% of the Matthieu equation we have Sys_dim=2;
Sys_dim=2;

% Constant defined in eq. (15)
% Exemplary calculation for the Matthieu equation
% d A /d phi_1=[0 0 ; eps*ampls(1)*sin(Om_1*t+phi_1) 0]
% d A /d phi_2=[0 0 ; eps*ampls(2)*sin(Om_2*t+phi_2) 0]
% => sup || d A /d phi_1+d A /d phi_1||=sup
% |eps*ampls(1)*sin(phi_1)+eps*ampls(2)*sin(phi_2)|<2*eps
CA=2*eps;

% initialize stability map
stab=zeros(length(ks),length(ws));


I=eye(Sys_dim);
% integrating equations of variation
% sweep parameter space

pars=[0.2 c eps.*ampls];
Oms=[1 sqrt(3)/3];

[t,x]= ode45(@(t,x)eq_of_var(t,x,Oms,[0;0],pars),[0 800],I(:));
figure
subplot(2,1,1)
plot(t,vecnorm(x(:,3:4).'))
del=0.1;
p1=patch([min(t)-del max(t)+del max(t)+del min(t)-del],...
     [min(vecnorm(x(:,3:4).'))-del min(vecnorm(x(:,3:4).'))-del max(vecnorm(x(:,3:4).'))+del max(vecnorm(x(:,3:4).'))+del],...
    'white');
set(p1,'Facealpha',0, 'Edgecolor','green','Linewidth',2 )
axis tight
title(['Parameters: $k=0.2$ and $\Omega=\sqrt{3}/3$  ' ],'Fontsize',22,'Interpreter','latex')
 xlabel('time','Fontsize',22,'Interpreter','latex')
ylabel(' $|x(t)|$','Fontsize',22,'Interpreter','latex')
set(gca,'fontsize',22)


pars=[0.21 c eps.*ampls];
Oms=[1 sqrt(3)/3];

[t,x]= ode45(@(t,x)eq_of_var(t,x,Oms,[0;0],pars),[0 800],I(:));
subplot(2,1,2)

plot(t,vecnorm(x(:,3:4).'))

 p1=patch([min(t)-del max(t)+del max(t)+del min(t)-del],...
     [min(vecnorm(x(:,3:4).'))-del min(vecnorm(x(:,3:4).'))-del max(vecnorm(x(:,3:4).'))+del max(vecnorm(x(:,3:4).'))+del],...
    'white');
set(p1,'Facealpha',0, 'Edgecolor','m','Linewidth',2 )
%set(gca,'XTick',
axis tight

title(['Parameters: $k=0.21$ and $\Omega=\sqrt{3}/3$  ' ],'Fontsize',22,'Interpreter','latex')
xlabel('time','Fontsize',22,'Interpreter','latex')
ylabel(' $|x(t)|$','Fontsize',22,'Interpreter','latex')
set(gca,'fontsize',22)
set(gcf,'Position',[ 500   55   800   500])

