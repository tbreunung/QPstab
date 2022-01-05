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
Nmax=50;
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

% sweep parameter space
for ii=1:          length(ks)
    ii
    
    
    for jj=1 :    length(ws)
        jj
        tic
          pars=[ks(ii) c eps.*ampls];
        Oms=[1 ws(jj)];
        
        
        % call algoritm to assess the stability of the origin of system (1) 
        tmp=cover_torus(pars,Sys_dim,Oms,CA,Nmax,delta_tol,N_edge,nw);
        %save result in the statbility map
        stab(ii,jj)=tmp;
        toc
    end
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
save('stab_map_matthieu')

%%
function stab=cover_torus(pars,Sys_dim,Oms,CA,Nmax,delta_tol,N_edge,nw)
% function carrying out the stability investigations of a quasi-periodic
% orbit
%
% input:    pars       parameter for equations of variations
%           Sys_dim    Dimension of system 
%           Oms        Frequency vector
%           CA         Constant bounding the matrix A (cf. eq. (15))
%           Nmax       Maximal number to periods for which we integrate the
%                      equations of variation
%                      i.e. Tmax=Nmax*T=Nmax*2*pi/Oms(1)
%           delta_tol  tolerance, if estimated neighbourhood, i.e. 
%                      Delta(phi_0) shrinks below delta_tol we terminate 
%                      our calculations
%           N_edge     estimated number of neighbourhoods along each edge
%                      of the cubic section of the trap
%                      MUST be a positive integer, 
%                      we set N_edge=10 in all our calculations
%           nw         number of workers
%                      depends on the machine used, MUST be a positive integer
%
% output: stab        -2: overflow in equation of variation
%                     -1: neighbourhood, i.e. deltla, below tolearance
%                      0: no contraction within maximal time span
%                         observed
%                      1: orbit is asymptotically stable

% creating identitiy matrix 
I=eye(Sys_dim);

% initializing stab
stab=1;
% dimension of the trap
Freq_dim=length(Oms)-1;

% sampling time T 
Ts=2*pi/Oms(1);


% We obtain the neighbourhood U(phi) for the initial anlge phi=0 to get an
% idea how many neighbourhoods are neccessary to cover the trap.

% integrate equations of variations form the initial angle  phi=0 until Tmax
[t,x]=ode45(@(t,x)eq_of_var(t,x,Oms,zeros(1,Freq_dim+1),pars),0:Ts:Ts*Nmax,I(:)); 

% initializing rho_max and rho_min
rho_max=zeros(length(t),1);
rho_min=zeros(length(t),1);
 
for ii=1:length(t)
    
    % Calculate PHI*PHI.'
    PHI=reshape(x(ii,:),Sys_dim,Sys_dim);
    CG=PHI*PHI.';

    % check for overflow 
    if max(isnan(CG(:)))==1 || max(isinf(CG(:)))==1  
        % display overflow result
        disp(['Overflow observed at angle phi=' num2str(zeros(1,Freq_dim+1))])
        % set stab to two
        stab=-2;
        % terminate calculations
        break
        
    else
        % calculate and assign rho_max and rho_min
        CG_eigs=eig(CG);
        rho_max(ii)=sqrt(max(CG_eigs));
        rho_min(ii)=sqrt(min(CG_eigs));
    end
end
% continue only if stability calculations are not terminated
if  stab==1
    % check if rho_max is less than one for some NT
    if  min(rho_max)<1
        % calculate eta
        eta=cumtrapz(t,rho_max./rho_min );
        
        % calculate neighbourhood and display it  
        delta_max=max(-log(rho_max)./(CA.*eta))
        % check is neighbourhood is less than tolerance specified
        if delta_max<delta_tol
            % if delta is less than tolearnce, display result and set stab
            % to negative one
            disp(['Radius of ball to small, caclulations terminated delta=' num2str(delta_max)])
            stab=-1;
        end
    else
        % if rho_max>1 for all N less than Nmax, display result and
        % set stab to zero
        disp(['No contraction observed at initial angle phi='  num2str(zeros(1,Freq_dim+1) )])
        stab=0;
    end
end
% proceed only if stability calculations are not terminated
if stab==1
    % We devided the trap into cubic sections with edge length
    % N_edge*delta_max
    N_sec=floor(2*pi/(delta_max*N_edge));
    % if N_sec<nw we divide the trap into nw(=number of workers) cubic
    % sections
    if N_sec<nw
        N_sec=nw;
    end
    % Dispaly an estimate of how many times the equations of variations
    % need to be integrated
    N_sec^Freq_dim
    % Initialize sections
    grid=linspace(0,2*pi,N_sec);
    % Sectioning the trap
    for ii=1:Freq_dim-1
        grid=combvec(grid,linspace(0,2*pi,N_sec));
    end
    % Edge length of each trap
    del_phi=2*pi/(N_sec-1);
    
    % Distribute sections to the workers 
    % Each workers has to cover about bx sections
    bx=floor(length(grid(1,:))/nw);
    % Remainder sections are distributed among the workers
    remainder=mod(length(grid(1,:)),nw);
    vec1=bx.*(1:nw);
    vec2=remainder.*ones(1,nw);
    vec2(1:remainder)=1:remainder;
    % idx_end is the index of the last section each individual worker 
    % has to cover 
    idx_end=vec1+vec2;
    % idx_start is the index of the first section each individual worker
    % has to cover
    idx_start=ones(1,nw)+idx_end;
    idx_start(end)=[];
    idx_start =[1 idx_start];
     
    % start parallel execution with cross communication 
    spmd
        % number of worker
        Par_idx=labindex;
        % Each worker covers the sections with index idx_start(Par_idx) 
        % until idx_end(Par_idx)
        for box_idx=idx_start(Par_idx):idx_end(Par_idx)
            % Exemplary call to monitor progress of stability 
            % investigations (usually deactivated to avoid exessive output)    
            %       disp([ 'Worker ' num2str(Par_idx) ' : ' num2str(floor(10^4*(box_idx-idx_start(Par_idx))/(idx_end(Par_idx)-idx_start(Par_idx)))/100) '% complete     stab=' num2str(stab) ])
          
            % Check if other worker has terminated calculations
            if stab==1
                if labProbe==1
                    % if some other worker send an update on stab variable
                    % this update is received
                    stab=labReceive;
                    
                end
            end
            % proceed only if stability calculations are not terminated 
            if stab==1
                
                % all initial angles are biased by the corresponding
                % section angle
                phi0=grid(:,box_idx).';
                % initialize current initial angle
                phi=zeros(1,Freq_dim);
                % initialize angle for incrementing
                phi_step=zeros(1,Freq_dim);
                while phi(end)<=del_phi
                    % call to integrate equation of variation 
                    % integration is started from the initial angle
                    % phi0+phi
                     
                    [t,x]=call_ode(pars,[0 phi0+phi],Oms,Nmax,Sys_dim);
                    
                    % initializing rho_max and rho_min
                    rho_max=zeros(length(t),1);
                    rho_min=rho_max;
 
                    for ii=1:length(t)
                        % Calculate PHI*PHI.'
                        PHI= reshape(x(ii,:),Sys_dim,Sys_dim);
                        CG=PHI*PHI.';
                        % check for overflow
                         if max(isnan(CG(:)))==1 || max(isinf(CG(:)))==1   
                            %disp(['Overflow observed at angle phi=' num2str([0 phi0+phi])])
                            % set stab to negative two if overflow is
                            % detected
                            stab=-2;
                            % send updated stab to the other workers
                            labSend(stab,[1:Par_idx-1 Par_idx+1:nw]);
                            % terminate calculations
                            break
                         else
                            % assign values to rho_max and rho_min
                            CG_eigs=eig(CG);
                            
                            rho_max(ii)=sqrt(max(CG_eigs));
                            rho_min(ii)=sqrt(min(CG_eigs));
                        end
                    end
                    % proceed only if stability calculations are not terminated 
                    if stab==1
                        % check if rho_max is less than one for some NT           
                        if   min(rho_max)<1
                            % calculate eta
                             eta=cumtrapz(t,rho_max./rho_min );
                            % calculate neighbourhood/ delta
                            delta_max=max(-log(rho_max)./(CA.*eta));

                   
                            % check if delta is less than specified
                            % tolerance
                            if delta_max<delta_tol
                                 %   disp(['Radius of ball to small, caclulations terminated delta=' num2str(delta_max )])
                                % set stab to negative one if neighbourhood
                                % is below tolerance
                                stab=-1;
                                % send updated stab to the other workers
                                labSend(stab,[1:Par_idx-1 Par_idx+1:nw]); %  
                                % terminate calculations
                                break
                            else
                                % if delta_max it greater than tolerance,
                                % we set the first entry of the the vector 
                                % phi_step to delta_max
                                % phi_step will be used for incrementing
                                % phi to cover the trap with neighbourhoods
                                phi_step(1)=delta_max;
                                % all zero entries of phi_step are set to
                                % delta_max
                                phi_step(phi_step(1:end)==0)=delta_max;
                                % non-zero entries of phi_step are compared
                                % to delta_max and phi_step is updated to 
                                % be the minimum value   
                                phi_step(phi_step(1:end)>0)=min(phi_step(1:end),delta_max.*ones(1,Freq_dim));
                                % check if phi+phi_step exceed the
                                % neighbourhood
                                kk=1;
                                while kk<Freq_dim
                                    % if phi(kk)+delta_max(kk) exceed delta_max 
                                    % (=edge length of trap section) then
                                    % the next phase is incremented
                                    if phi(kk)+delta_max>del_phi
                                        kk=kk+1;
                                    else
                                        break
                                    end
                                    
                                end
                                % the kk-th phase is incremented by
                                % phi_step(kk) 
                                phi(kk)=phi(kk)+phi_step(kk)+delta_tol;
                                % phases below kk are set to zero to
                                % restart the incrementing process
                                phi(1:kk-1)=0.*phi(1:kk-1);
                                phi_step(1:kk-1)=0.*phi_step(1:kk-1);
                            end
                            
                            
                        else
                            % if rho_max is greater than one for all Nmax,
                            % then stab is set to zero
                            stab=0;
                            % send updated stab to other workers
                            labSend(stab,[1:Par_idx-1 Par_idx+1:nw]);
                            % terminate calculations
                            break
                        end
                    else
                        % termniate calculations of stab is not equal to
                        % one (i.e. stability investigations are
                        % inconclusive)
                        break
                    end
                end
                
                
            end
        end
       
    end
    % set output to be the minimum of all stab variables from all workers
    stab=min([stab{:} ]);
end

end


%%

function [t,x]=call_ode(pars,phi,Oms,Nmax,Sys_dim)
% Call to integrate equation of variation

% Identity matrix as initial condition
I=eye(Sys_dim);
% sampling time Ts
Ts=2*pi/Oms(1);
% integrating equations of variation
[t,x]= ode45(@(t,x)eq_of_var(t,x,Oms,phi,pars),0:Ts: Ts*Nmax,I(:));
 
end
 