% Exemplary call of QPstab to investigate the stability of the origin of
% the Matthieu equation

% q_dd+c*q_d+[k+eps*(ampls(1)*cos(Om_1 *t)+ampls(1)*cos(Om_2 *t) ) ] q=0


% First order form of the equation of variation must be as separate function
% in a m-file (e.g. eq_of_var.m for the above example)
clear all
close all


% Code to run QPstab on ETHs Euler cluster (acess restricetd to ETH
% affiliates)
%mode='Euler';
%if strcmp(mode,'Euler')==1
%   batch_job = parcluster('EulerLSF8h');
%    nw=144;
    
%    batch_job.SubmitArguments = '-W 24:00 -R "rusage[mem=2048]"';
%    pool = parpool(batch_job,nw);
%else
%    nw=16;
    

% start parallel computing pool if not already started  
if isempty(gcp('nocreate'))==1
    pool=parpool;
    % get number of workers
    nw = pool.NumWorkers;
else 
    tmp=gcp;
    nw=tmp.NumWorkers;
end




% We will vary Om_2 and the stiffness k
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
for ii=1:          length(ks)
    ii
    
     tic
    for jj=1 :    length(ws)
         
        pars=[ks(ii) c eps.*ampls];
        Oms=[1 ws(jj)];
        Ts=2*pi/Oms(1);
        
        [t,x]= ode45(@(t,x)eq_of_var(t,x,Oms,[0;0],pars),0:Ts: Ts*Nmax,I(:));
        PHI= reshape(x(end,:),Sys_dim,Sys_dim);
        CG=PHI*PHI.';
        FTLE=1/(2*Ts*Nmax)*log(max(eig(CG)));
        
        if FTLE<0
            stab(ii,jj)=1;
        end
        
        
        
    end
    toc
    % exemplary commad to save intermediate results
  %  save(['stab_map_matthieu_inter_ii=' num2str(ii)],'ws','ks','eps','ampls','Nmax','N','delta_tol','Sys_dim','CA','stab','ii','jj')
    
end

        

%
% if strcmp(mode,'Euler')==1
%     pool.delete()
%     clear batch_job
% end
%


% save results 
save('stab_map_matthieu_FTLE')

%%
 