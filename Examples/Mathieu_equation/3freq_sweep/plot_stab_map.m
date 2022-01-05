clear all
close all
%load('stab_map_matthieu3D.mat')

load('stab_map_matthieu3D_eps02_15^3.mat')

stab_map=figure;
hold on
%stab(ii+1:end)=[];
%set(stab_map,'Visible','off')
idx=find(stab==1);
if isempty(idx)==0
pl(1,:)=plot3(OMS_grid(1,idx),OMS_grid(2,idx),OMS_grid(3,idx),'dg','MarkerFaceColor','g');
        plot3(OMS_grid(1,idx),OMS_grid(3,idx),OMS_grid(2,idx),'dg','MarkerFaceColor','g');
        plot3(OMS_grid(2,idx),OMS_grid(1,idx),OMS_grid(3,idx),'dg','MarkerFaceColor','g');
        plot3(OMS_grid(2,idx),OMS_grid(3,idx),OMS_grid(1,idx),'dg','MarkerFaceColor','g');
        plot3(OMS_grid(3,idx),OMS_grid(1,idx),OMS_grid(2,idx),'dg','MarkerFaceColor','g');
        plot3(OMS_grid(3,idx),OMS_grid(2,idx),OMS_grid(1,idx),'dg','MarkerFaceColor','g');
end

idx=find(stab==0);
if isempty(idx)==0
pl(2,:)=plot3(OMS_grid(1,idx),OMS_grid(2,idx),OMS_grid(3,idx),'om','MarkerFaceColor','m');
        plot3(OMS_grid(1,idx),OMS_grid(3,idx),OMS_grid(2,idx),'om','MarkerFaceColor','m');
        plot3(OMS_grid(2,idx),OMS_grid(1,idx),OMS_grid(3,idx),'om','MarkerFaceColor','m');
        plot3(OMS_grid(2,idx),OMS_grid(3,idx),OMS_grid(1,idx),'om','MarkerFaceColor','m');
        plot3(OMS_grid(3,idx),OMS_grid(1,idx),OMS_grid(2,idx),'om','MarkerFaceColor','m');
        plot3(OMS_grid(3,idx),OMS_grid(2,idx),OMS_grid(1,idx),'om','MarkerFaceColor','m');

end

idx=find(stab==-1);
if isempty(idx)==0
pl(3,:)=plot3(OMS_grid(1,idx),OMS_grid(2,idx),OMS_grid(3,idx),'sk','MarkerFaceColor','k');
        plot3(OMS_grid(1,idx),OMS_grid(3,idx),OMS_grid(2,idx),'sk','MarkerFaceColor','k');
        plot3(OMS_grid(2,idx),OMS_grid(1,idx),OMS_grid(3,idx),'sk','MarkerFaceColor','k');
        plot3(OMS_grid(2,idx),OMS_grid(3,idx),OMS_grid(1,idx),'sk','MarkerFaceColor','k');
        plot3(OMS_grid(3,idx),OMS_grid(1,idx),OMS_grid(2,idx),'sk','MarkerFaceColor','k');
        plot3(OMS_grid(3,idx),OMS_grid(2,idx),OMS_grid(1,idx),'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
pl(4,:)=plot3(OMS_grid(1,idx),OMS_grid(2,idx),OMS_grid(3,idx),'xr');
        plot3(OMS_grid(1,idx),OMS_grid(3,idx),OMS_grid(2,idx),'xr');
        plot3(OMS_grid(2,idx),OMS_grid(1,idx),OMS_grid(3,idx),'xr');
        plot3(OMS_grid(2,idx),OMS_grid(3,idx),OMS_grid(1,idx),'xr');
        plot3(OMS_grid(3,idx),OMS_grid(1,idx),OMS_grid(2,idx),'xr');
        plot3(OMS_grid(3,idx),OMS_grid(2,idx),OMS_grid(1,idx),'xr');
end
 
xlabel('$\Omega_1$','Fontsize',16,'Interpreter','latex')
ylabel('$\Omega_2$','Fontsize',16,'Interpreter','latex')
zlabel('$\Omega_3$','Fontsize',16,'Interpreter','latex')

%leg=legend(pl,[  'stable' ],['~~' newline '~~' newline 'no contraction' newline 'observed' newline '~~'],['$\Delta$ below'  newline 'tolerance']);
%    '$\Omega_j=2$','$\Omega_1+\Omega_2=2 $  ','$\Omega_1-\Omega_2=2$ ','$2\Omega_1+\Omega_2=2$',...
%    '$-2\Omega_1+ \Omega_2=2$ ','$\Omega_1-2 \Omega_2=2$ ');
%set(leg,'Fontsize',16,'location','NorthEastOutside','Interpreter','latex')
%set(leg,'Box','off')
del=0.3;
p1=patch([min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del max(OMS_grid(1,:))+del min(OMS_grid(1,:))-del],...
    sqrt(5)/2.*[min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del max(OMS_grid(1,:))+del min(OMS_grid(1,:))-del],...
    [min(OMS_grid(3,:))-del min(OMS_grid(3,:))-del max(OMS_grid(3,:))+del max(OMS_grid(3,:))+del],...
    'blue')
set(p1,'Facealpha',0.5,'Edgecolor','blue','Linewidth',3 )
p2=patch([min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del max(OMS_grid(1,:))+del min(OMS_grid(1,:))-del],...
     [min(OMS_grid(2,:))-del min(OMS_grid(2,:))-del max(OMS_grid(2,:))+del max(OMS_grid(2,:))+del],...
    sqrt(3).*[1 1 1 1],...
    'red')
set(p2,'Facealpha',0.5 ,'Edgecolor','red','Linewidth',3 )


set(gca,'fontsize',16)
set(gcf,'Position',[ 500   55   500   500])
set(gca,'View',[ -24.8221   15.5281])
del=del*1.1;
axis([min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del min(OMS_grid(2,:))-del max(OMS_grid(2,:))+del min(OMS_grid(3,:))-del max(OMS_grid(3,:))+del])
%set(stab_map,'Visible','on')
%xlabel('\Omega_1')
%ylabel('\Omega_2')
%zlabel('\Omega_3')
%leg=legend(pl,'stable','no contraction observed' ,'terminated \delta below tol');
%set(leg,'location','NorthEastOutside')