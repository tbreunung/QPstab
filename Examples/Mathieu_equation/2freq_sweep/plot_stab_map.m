close all
clear all
load('stab_map_matthieu100x100.mat')

stab_map=figure;
hold on

sz=2;

idx=find(stab==0);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(2,:)=plot(w1(ridx),w2(cidx),'om','MarkerFaceColor','m','MarkerSize',sz);
pl(2,:)=plot(100,100,'om','MarkerFaceColor','m');

end

idx=find(stab==-1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(3,:)=plot(w1(ridx),w2(cidx),'sk','MarkerFaceColor','k','MarkerSize',sz);
pl(3,:)=plot(100,100,'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(4,:)=plot(w1(ridx),w2(cidx),'xr');
end

idx=find(stab==1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(1,:)=plot(w1(ridx),w2(cidx),'dg','MarkerFaceColor','g','MarkerSize',sz);
pl(1,:)=plot(100,100,'dg','MarkerFaceColor','g');

end
Nx=length(w1);
pl(4,:)=plot([2 2],[0 3],'--k');
plot([0:0.5: 3],[2 2 2 2 2 2 2],'--k');

plot(w1,(2.*ones(1,Nx)-w1)./1,'-k')
%plot((2.*ones(1,100)-w2)./1,w2,'-k')
pl(5,:)=plot(w1(1:5:end),(2.*ones(1,length(w1(1:5:end)))-w1(1:5:end))./1,'-dk');
%plot((2.*ones(1,length(w1(1:5:100)))-w2(1:5:100))./1,w2(1:5:100),'-dk')

plot(w1,(2.*ones(1,Nx)+w1)./1,'-k');
plot((2.*ones(1,Nx)+w2)./1,w2,'-k')
pl(6,:)=plot(w1(1:5:end),(2.*ones(1,length(w1(1:5:end)))+w1(1:5:end))./1,'-sk');
plot((2.*ones(1,length(w2(1:5:end)))+w2(1:5:end))./1,w2(1:5:end),'-sk')

plot(w1,(2.*ones(1,Nx)-2.*w1)./1,'-k');
plot((2.*ones(1,Nx)-2.*w2)./1,w2,'-k')
pl(7,:)=plot(w1(1:5:end),(2.*ones(1,length(w1(1:5:end)))-2.*w1(1:5:end))./1,'-ok');
plot((2.*ones(1,length(w2(1:5:end)))-2.*w2(1:5:end))./1,w2(1:5:end),'-ok');


pl(8,:)=plot(w1,(2.*ones(1,Nx)+2.*w1)./1,'-.k');


plot((2.*ones(1,Nx)+2.*w2)./1,w2,'-k');
pl(9,:)=plot((2.*ones(1,length(w2(1:5:end)))+2.*w2(1:5:end))./1,w2(1:5:end),'-^k');

 

xlabel('$\Omega_1$','Fontsize',16,'Interpreter','latex')
ylabel('$\Omega_2$','Fontsize',16,'Interpreter','latex')
axis([w1(1) w1(end) w2(1) w2(end)])
leg=legend(pl,'asymptotically stable origin','no contraction observed' ,'$\Delta$ below tolerance',...
    '~~$\Omega_j=2$','~~$\Omega_1+\Omega_2=2 $  ','~~$\Omega_1-\Omega_2=2$ ','~~$2\Omega_1+\Omega_2=2$',...
    '~~$-2\Omega_1+ \Omega_2=2$ ','~~$\Omega_1-2 \Omega_2=2$ ');
set(leg,'Fontsize',16,'location','NorthEastOutside','Interpreter','latex')
set(leg,'Box','off')
set(gca,'fontsize',16)
set(gcf,'Position',[ 500   55   800   500])