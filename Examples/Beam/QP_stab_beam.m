
clear all
close all

%parallel.importProfile('/cluster/apps/matlab/support/EulerLSF8h.settings')


mode='Euler';
if strcmp(mode,'Euler')==1
    batch_job = parcluster('EulerLSF8h');
    nw=144;
    %delete(batch_job.Jobs);
    batch_job.SubmitArguments = '-W 24:00 -R "rusage[mem=2048]"';
    pool = parpool(batch_job,nw);
else
    nw=16;
    
    if isempty(gcp('nocreate'))==1
        pool=parpool;
        
        %c = parcluster('local');
        nw = pool.NumWorkers;
        
    end
end

Ndis=10;
Sys_dim=Ndis-4;
Freq_dim=2;

%L=10;


Nmax=100;
N=10;
delta_tol= 1*10^-4;

r=0.1;
mu=0.01;

Oms=[1 sqrt(3)];

C=2*mu.*eye(Sys_dim);


% OM_low=0.1;
% OM_high=3;
% Nfreq=50;
% 
% OMS_grid=linspace(OM_low,OM_high,Nfreq);
% for ii=1:Freq_dim-1
%     OMS_grid=combvec(OMS_grid,linspace(OM_low,OM_high,Nfreq));
%     [~,idx] = unique(sort(OMS_grid.',2),'rows','stable');
%     OMS_grid = OMS_grid(:,idx);
% end
% irr=10^-5.*sqrt(nthprime(1:Freq_dim)).'*ones(1,length(OMS_grid(1,:)));
% OMS_grid=OMS_grid+irr;

Npars=50;
L=linspace(8,14,Npars);
eps=linspace(1,7.5,Npars).*10^-3;

%OM_3=linspace(OM_low,OM_high,Nfreq);
%%
stab=zeros(Npars,Npars);

for ii=1:Npars
    ii
    for jj=1:Npars
    jj
    h=L(ii)/Ndis;
    K_2=1/(r^2*h^2).*spdiags(ones(Sys_dim,1)*[1 -2 1],-1:1,Sys_dim,Sys_dim);
    K_1=1/h^4.*spdiags(ones(Sys_dim,1)*[1 -4 6 -4 1],-2:2,Sys_dim,Sys_dim);
     
    ampls=eps(jj).*ones(1,Freq_dim);
    
    CA=max(abs(eig(K_2)))*Freq_dim*eps(jj);
    %Oms=OMS_grid(:,ii);
    
    tic
    
    %tmp=cover_torus(@(t,x,phi)eq_of_var(t,x,Oms,phi,pars),Tmax,N,CA,delta_tol,Sys_dim,Freq_dim,Oms,eps,nw);
    tmp=cover_torus(C,K_1,K_2,ampls,Oms,CA,Nmax,delta_tol,N,nw);
    stab(ii,jj)=tmp;
    toc
    save('stat_beam','ii','jj', 'L','eps' ,'stab')
    end
end




if strcmp(mode,'Euler')==1
    pool.delete()
    clear batch_job
    save('stab_map_beam')
end


%%
function stab=cover_torus(C,K_1,K_2,ampls,Oms,CA,Nmax,delta_tol,N,nw)
Sys_dim=2*min(size(C));
I=eye(Sys_dim);
stab=1;

Freq_dim=length(Oms)-1;

Ts=2*pi/Oms(1);
[t,x]=ode45(@(t,x)Beam_eq_of_var(t,x,Oms,ampls,C,K_1,K_2,zeros(1,Freq_dim+1)),0:Ts:Ts*Nmax,I(:));%,odeset(  'RelTol',1e-10,'AbsTol',1e-10));

eig_max=zeros(length(t),1);
eig_min=zeros(length(t),1);

for ii=1:length(t)
    
    PHI=reshape(x(ii,:),Sys_dim,Sys_dim);
    CG=PHI*PHI.';
    
    
    if max(isnan(CG(:)))==1 || max(isinf(CG(:)))==1
        disp(['Overflow observed at angle phi=' num2str(zeros(1,Freq_dim+1))])
        stab=-2;
        break
    else
        CG_eigs=eig(CG);
        eig_max(ii)=sqrt(max(CG_eigs));
        tmp=min(CG_eigs);
        if tmp< 0 %tol_eig
            eig_min(ii:end)=[];
            eig_max(ii:end)=[];
            t(ii:end,:)=[];
            break
        else
            eig_min(ii)=sqrt(tmp);
        end
        
    end
end
if  stab==1
    
    if  min(eig_max)<1
        rho_tilde=cumtrapz(t,eig_max./eig_min );
        
        %CAA=ampl.*eig_max.*cumtrapz(t,1./eig_min );
        delta_max=max(-log(eig_max)./(CA.*rho_tilde))
        
        %delta_max=max(-log(eig_max)./C2)
        %delta_max= (log(Cphi_ini)-alpha_ini*t(end))/(Cphi_ini*CAA)
        %         [val,idx]=min(eig_max);
        %
        if delta_max<delta_tol
            disp(['Radius of ball to small, caclulations terminated delta=' num2str(delta_max)])
            stab=-1;
            %              figure
            %         plot(t,eig_max)
        end
    else
        disp(['No contraction observed at initial angle phi='  num2str(zeros(1,Freq_dim+1) )])
        %                 figure
        %                plot(t,eig_max)
        stab=0;
    end
end
if stab==1
    N_sec=floor(2*pi/(delta_max*N));
    if N_sec<nw
        N_sec=nw;
    end
    N_sec^Freq_dim
    grid=linspace(0,2*pi,N_sec);
    for ii=1:Freq_dim-1
        grid=combvec(grid,linspace(0,2*pi,N_sec));
    end
    %[phi1, phi2]=ndgrid(linspace(0,2*pi,N*nw));
    del_phi=2*pi/(N_sec-1);
    
    bx=floor(length(grid(1,:))/nw);
    remainder=mod(length(grid(1,:)),nw);
    vec1=bx.*(1:nw);
    vec2=remainder.*ones(1,nw);
    vec2(1:remainder)=1:remainder;
    idx_end=vec1+vec2;
    idx_start=ones(1,nw)+idx_end;
    idx_start(end)=[];
    idx_start =[1 idx_start];
    %
    
    spmd
        Par_idx=labindex;
        for box_idx=idx_start(Par_idx):idx_end(Par_idx)
            %       disp([ 'Worker ' num2str(Par_idx) ' : ' num2str(floor(10^4*(box_idx-idx_start(Par_idx))/(idx_end(Par_idx)-idx_start(Par_idx)))/100) '% complete     stab=' num2str(stab) ])
            if stab==1
                if labProbe==1
                    
                    stab=labReceive;
                    
                end
            end
            
            if stab==1
                
                
                phi0=grid(:,box_idx).';
                phi=zeros(1,Freq_dim);
                phi_step=zeros(1,Freq_dim);
                while phi(end)<=del_phi
                    %   phi0+phi
                    [t,x]=call_ode(C,K_1,K_2,ampls,[0 phi0+phi],Oms,Nmax);
                    
                    eig_max=zeros(length(t),1);
                    eig_min=eig_max;
                    %  eig_max2=eig_max;
                    
                    % Phi_t0t1= reshape(x(end,:),Sys_dim,Sys_dim);
                    
                    for ii=1:length(t)
                        %  PHI2=Phi_t0t1*inv(reshape(x(ii,:),Sys_dim,Sys_dim));  %reshape(x(ii,:),Sys_dim,Sys_dim);%      %Phi_t0t1*inv([x(ii,1:2).'  x(ii,3:4).']);%%
                        PHI= reshape(x(ii,:),Sys_dim,Sys_dim);
                        CG=PHI*PHI.';
                        %  CG2=PHI2*PHI2.';
                        
                        if max(isnan(CG(:)))==1 || max(isinf(CG(:)))==1
                            disp(['Overflow observed at angle phi=' num2str([0 phi0+phi])])
                            stab=-2;
                            labSend(stab,[1:Par_idx-1 Par_idx+1:nw]);
                            break
                        else
                            CG_eigs=eig(CG);
                            
                            eig_max(ii)=sqrt(max(CG_eigs));
                            tmp=min(CG_eigs);
                            if tmp< 0 %tol_eig
                                eig_min(ii:end)=[];
                                eig_max(ii:end)=[];
                                t(ii:end,:)=[];
                                break
                            else
                                eig_min(ii)=sqrt(tmp);
                            end
                            
                        end
                    end
                    if stab==1
                        
                        if   min(eig_max)<1
                            rho_tilde=cumtrapz(t,eig_max./eig_min );
                            
                            delta_max=max(-log(eig_max)./(CA.*rho_tilde));
                            
                            %        CAA=ampl.*eig_max.*cumtrapz(t,1./eig_min);
                            %       delta_max=max(-log(eig_max)./CAA);
                            
                            
                            if delta_max<delta_tol
                                %   disp(['Radius of ball to small, caclulations terminated delta=' num2str(delta_max )])
                                
                                stab=-1;
                                %if labProbe==0
                                labSend(stab,[1:Par_idx-1 Par_idx+1:nw]); % Par_idx+1:nw gop(@min,stab,'all')
                                % end
                                % stab=-1;
                                break
                            else
                                phi_step(1)=delta_max;
                                phi_step(phi_step(1:end)==0)=delta_max;
                                phi_step(phi_step(1:end)>0)=min(phi_step(1:end),delta_max.*ones(1,Freq_dim));
                                kk=1;
                                while kk<Freq_dim
                                    if phi(kk)+delta_max>del_phi
                                        kk=kk+1;
                                    else
                                        break
                                    end
                                    
                                end
                                
                                phi(kk)=phi(kk)+phi_step(kk)+delta_tol;%*(sqrt(3)/2+eq(kk,1)*(1-sqrt(3)/2));
                                
                                phi(1:kk-1)=0.*phi(1:kk-1);
                                phi_step(1:kk-1)=0.*phi_step(1:kk-1);
                            end
                            
                            
                        else
                            %   disp(['No contraction observed at phi='  num2str((phi0+phi) )])
                            stab=0;
                            %if labProbe==0
                            labSend(stab,[1:Par_idx-1 Par_idx+1:nw]);
                            % end
                            %stab=0;
                            break
                        end
                    else
                        break
                    end
                end
                
                
            end
        end
        
    end
    stab=min([stab{:} ]);
end

end


%set(f,'Visible','on')
%%

function [t,x]=call_ode(C,K_1,K_2,ampls,phi,Oms,Nmax)
I=eye(2*min(size(C)));
Ts=2*pi/Oms(1);

[t,x]= ode45(@(t,x)Beam_eq_of_var(t,x,Oms,ampls,C,K_1,K_2,phi),0:Ts: Ts*Nmax,I(:));
%[t,x]= ode45(@(t,x)my_eq_var(t,x,phi),[0 Tmax],I(:));
end

