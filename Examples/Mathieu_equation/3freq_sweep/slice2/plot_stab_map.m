clear all
%close all
load('stab_map_matthieu3D_slice2.mat')
 
stab_map=figure;
hold on



%set(stab_map,'Visible','off')

idx=find(stab==0);
if isempty(idx)==0
pl(2,:)=plot(OMS_grid(1,idx),OMS_grid(2,idx),'om','MarkerFaceColor','m'); 
plot(OMS_grid(2,idx),OMS_grid(1,idx),'om','MarkerFaceColor','m'); 
end

idx=find(stab==-1);
if isempty(idx)==0
pl(3,:)=plot(OMS_grid(1,idx),OMS_grid(2,idx),'sk','MarkerFaceColor','k');
 plot(OMS_grid(2,idx),OMS_grid(1,idx),'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
pl(4,:)=plot(OMS_grid(1,idx),OMS_grid(3,idx),'xr');
end

idx=find(stab==1);
if isempty(idx)==0
pl(1,:)=plot(OMS_grid(1,idx),OMS_grid(2,idx),'dg','MarkerFaceColor','g');
plot(OMS_grid(2,idx),OMS_grid(1,idx),'dg','MarkerFaceColor','g');
end


% for ii=1:length(ks)
%     
%     for jj=1:length(ws)
%         figure(stab_map)
%         switch stab(ii,jj)
%             
%             case 1
%                 pl(1,:)=plot(ks(ii),ws(jj),'dg','MarkerFaceColor','g');
%             case 0
%                 pl(2,:)=plot(ks(ii),ws(jj),'om','MarkerFaceColor','m');
%             case -1
%                 pl(3,:)=plot(ks(ii),ws(jj),'sk','MarkerFaceColor','k');
%             case -2
%                 pl(4,:)=plot(ks(ii),ws(jj),'xr');
%         end
%      end
% end
del=0.1;
p1=patch([min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del max(OMS_grid(1,:))+del min(OMS_grid(1,:))-del],...
     [min(OMS_grid(2,:))-del min(OMS_grid(2,:))-del max(OMS_grid(2,:))+del max(OMS_grid(2,:))+del],...
    'white');
set(p1,'Facealpha',0, 'Edgecolor','red','Linewidth',2 )

%set(stab_map,'Visible','on')
xlabel('$\Omega_1$','Fontsize',16,'Interpreter','latex')
ylabel('$\Omega_2$','Fontsize',16,'Interpreter','latex')
title('$\Omega_3 =\sqrt{3}$','Interpreter','latex','Fontsize',16)
leg=legend(pl,[ 'asymptotically' newline 'stable origin' ],['~~' newline '~~' newline 'no contraction' newline 'observed' newline '~~'],['$\Delta$ below'  newline 'tolerance']);
set(leg,'Fontsize',16,'location','NorthEastOutside','Interpreter','latex')
set(leg,'Box','off')
  
%leg=legend(pl,'stable','no contraction observed' ,'terminated \Delta below tol');
%set(leg,'location','NorthEastOutside')
set(gca,'Xtick',0:1:3,'Ytick',0:1:3)
set(gca,'fontsize',16)
set(gcf,'Position',[ 500   55   700   500])
axis([min(OMS_grid(1,:))-del max(OMS_grid(1,:))+del min(OMS_grid(2,:))-del max(OMS_grid(2,:))+del])