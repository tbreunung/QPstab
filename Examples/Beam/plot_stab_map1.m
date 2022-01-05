% clear all
load('stab_map_beam_50x50.mat')

stab_map=figure;
hold on

%set(stab_map,'Visible','off')
idx=find(stab==1);
sz=3;
if isempty(idx)==0
plot(OMS_grid(1,idx),OMS_grid(2,idx),'dg','MarkerFaceColor','g','Markersize',sz);
plot(OMS_grid(2,idx),OMS_grid(1,idx),'dg','MarkerFaceColor','g','Markersize',sz);
pl(1,:)=plot(100,100,'dg','MarkerFaceColor','g');
end

idx=find(stab==0);
if isempty(idx)==0
plot(OMS_grid(1,idx),OMS_grid(2,idx),'om','MarkerFaceColor','m','Markersize',sz);
plot(OMS_grid(2,idx),OMS_grid(1,idx),'om','MarkerFaceColor','m','Markersize',sz);
pl(2,:)=plot(100,100,'om','MarkerFaceColor','m');

end

idx=find(stab==-1);
if isempty(idx)==0
plot(OMS_grid(1,idx),OMS_grid(2,idx),'sk','MarkerFaceColor','k','Markersize',sz);
plot(OMS_grid(2,idx),OMS_grid(1,idx),'sk','MarkerFaceColor','k','Markersize',sz);
pl(3,:)=plot(100,100,'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(4,:)=plot(ks(ridx),ws(cidx),'xr');
end
 

%set(stab_map,'Visible','on')
xlabel('$\Omega_1$','Fontsize',22,'Interpreter','latex')
ylabel('$\Omega_2$','Fontsize',22,'Interpreter','latex')
axis([min(OMS_grid(:)) max(OMS_grid(:)) min(OMS_grid(:)) max(OMS_grid(:))])
leg=legend(pl,[ 'asymptotically' newline 'stable origin' ],['~~' newline '~~' newline 'no contraction' newline 'observed' newline '~~'],['$\Delta$ below'  newline 'tolerance']);
set(leg,'Fontsize',22,'location','NorthEastOutside','Interpreter','latex')
set(leg,'Box','off')
set(gca,'fontsize',22)
set(gcf,'Position',[ 500   55   800   500])

