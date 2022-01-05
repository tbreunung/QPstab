clear all
close all
%load('stab_map_matthieu100x100.mat')
load('stab_map_matthieu_FTLE.mat')

stab_map=figure;
hold on

%set(stab_map,'Visible','off')
sz=2;
idx=find(stab==0);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(2,:)=plot(ks(ridx),ws(cidx),'om','MarkerFaceColor','m','MarkerSize',sz);
pl(2,:)=plot(100,100,'om','MarkerFaceColor','m');
end

idx=find(stab==-1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(3,:)=plot(ks(ridx),ws(cidx),'sk','MarkerFaceColor','k','MarkerSize',sz);
pl(3,:)=plot(100,100,'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(4,:)=plot(ks(ridx),ws(cidx),'xr');
end

idx=find(stab==1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(1,:)=plot(ks(ridx),ws(cidx),'dg','MarkerFaceColor','g','MarkerSize',sz);
pl(1,:)=plot(100,100,'dg','MarkerFaceColor','g');

end



 
xlabel('$k$','Fontsize',22,'Interpreter','latex')
ylabel('$\Omega_1$','Fontsize',22,'Interpreter','latex')
axis([ks(1) ks(end) ws(1) ws(end)])
leg=legend(pl,[  'FTLE$<0$' ],['~~' newline '~~' newline 'FTLE$>0$'  newline '~~' newline]);
set(leg,'Fontsize',22,'location','NorthEastOutside','Interpreter','latex')
set(leg,'Box','off')
set(gca,'fontsize',22)
set(gcf,'Position',[ 500   55   800   500])
% set(leg,'location','NorthEastOutside')
% leg.Fontsize=18;
% leg.Interpreter='latex';
% leg.location='NorthEastOutside';
%