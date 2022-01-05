% clear all
% close all
 %load('stab_map_beam.mat')
load('stab_map_beam2.mat')
stab_map=figure;
hold on

%set(stab_map,'Visible','off')
sz=3;
idx=find(stab==0);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(2,:)=plot(L(ridx),eps(cidx),'om','MarkerFaceColor','m','MarkerSize',sz);
pl(2,:)=plot(100,100,'om','MarkerFaceColor','m');
end

idx=find(stab==-1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(3,:)=plot(L(ridx),eps(cidx),'sk','MarkerFaceColor','k','MarkerSize',sz);
pl(3,:)=plot(100,100,'sk','MarkerFaceColor','k');

end

idx=find(stab==-2);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(4,:)=plot(L(ridx),eps(cidx),'xr');
end

idx=find(stab==1);
if isempty(idx)==0
[ridx,cidx]=ind2sub(size(stab),idx);
pl(1,:)=plot(L(ridx),eps(cidx),'dg','MarkerFaceColor','g','MarkerSize',sz);

end
pl(1,:)=plot(100,100,'dg','MarkerFaceColor','g');


% ws=zeros(length(L),Sys_dim);
% for ii=1:length(L)
% h=L(ii)/Ndis;
% K_1=1/h^4.*spdiags(ones(Sys_dim,1)*[1 -4 6 -4 1],-2:2,Sys_dim,Sys_dim);
% ws(ii,:)=sqrt(eigs(K_1,Sys_dim));
% end
%  plot(L,2.*ws,'-k');

xlabel(' Length $L$','Fontsize',22,'Interpreter','latex')
ylabel('Amplitude $a$','Fontsize',22,'Interpreter','latex')
axis([L(1) L(end) eps(1) eps(end)])
leg=legend(pl,[ 'asymptotically' newline 'stable origin' ],['~~' newline '~~' newline 'no contraction' newline 'observed' newline '~~'],['$\Delta$ below'  newline 'tolerance']);
 %   '$k=1/4$','$k=\frac{(\Omega_1+1)^2}{4} $  ','$k=\frac{(\Omega_1-1)^2}{4}$ ','$k=\frac{(2\Omega_1+1)^2}{4}$',...
 %   '$k=\frac{(2\Omega_1-1)^2}{4}$ ','$k=\frac{\Omega_1^2}{4}$ ','$k=~~\frac{ \Omega_1^2}{4}$');
 

set(leg,'Fontsize',22,'location','NorthEastOutside','Interpreter','latex')

set(leg,'Box','off')
% HeightScaleFactor = 1.5;
% NewHeight = leg.Position(4) * HeightScaleFactor;
% leg.Position(2) = leg.Position(2) - (NewHeight - leg.Position(4));
% leg.Position(4) = NewHeight;


set(gca,'fontsize',22)
set(gcf,'Position',[ 500   55   800   500])
% set(leg,'location','NorthEastOutside')
% leg.Fontsize=18;
% leg.Interpreter='latex';
% leg.location='NorthEastOutside';
%