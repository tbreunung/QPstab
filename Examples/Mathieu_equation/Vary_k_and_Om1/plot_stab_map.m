clear all
close all
load('stab_map_matthieu100x100.mat')
%load('stab_map_matthieu_FTLE.mat')

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


pl(4,:)=plot([0.25 0.25],[0 1.5],'--k');

plot((ws+ones(1,length(ks))).^2./4,ws,'-k');
pl(5,:)=plot((ws(1:2:100)+ones(1,length(ks(1:2:100)))).^2./4,ws(1:2:100),'-dk');


 plot((ws-ones(1,length(ks))).^2./4,ws,'- k');
pl(6,:)=plot((ws(1:5:100)-ones(1,length(ks(1:5:100)))).^2./4,ws(1:5:100),'-sk');

pl(7,:)=plot((2*ws+ones(1,length(ks))).^2./4,ws,'-ok');

pl(8,:)=plot((2*ws-ones(1,length(ks))).^2./4,ws,'-.k');

plot((ws+0.*ones(1,length(ks))).^2./4,ws,'-k');
pl(9,:)=plot((ws(1:5:100)+0.*ones(1,length(ks(1:5:100)))).^2./4,ws(1:5:100),'-xk');

plot((2*ws+0.*ones(1,length(ks))).^2./4,ws,'-k');
pl(10,:)=plot((2*ws(1:5:100)+0.*ones(1,length(ks(1:5:100)))).^2./4,ws(1:5:100),'-^k');
 
xlabel('$k$','Fontsize',16,'Interpreter','latex')
ylabel('$\Omega_1$','Fontsize',16,'Interpreter','latex')
axis([ks(1) ks(end) ws(1) ws(end)])
leg=legend(pl,'asymptotically stable origin','no contraction observed' ,'$\Delta$ below tolerance',...
    '~~$k=1/4$','$~~k=\frac{(\Omega_1+1)^2}{4} $  ','~~$k=\frac{(\Omega_1-1)^2}{4}$ ','~~$k=\frac{(2\Omega_1+1)^2}{4}$',...
    '~~$k=\frac{(2\Omega_1-1)^2}{4}$ ','~~$k=\frac{\Omega_1^2}{4}$ ','~~$k=~~\frac{ \Omega_1^2}{4}$');
set(leg,'Fontsize',16,'location','NorthEastOutside','Interpreter','latex')
set(leg,'Box','off')
set(gca,'fontsize',16)
set(gcf,'Position',[ 500   55   800   500])
% set(leg,'location','NorthEastOutside')
% leg.Fontsize=18;
% leg.Interpreter='latex';
% leg.location='NorthEastOutside';
%