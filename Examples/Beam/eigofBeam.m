L=linspace(8,14);
Ndis=10;
Sys_dim=Ndis-4;

for ii=1:length(L)
    h=L(ii)/Ndis;
    K_1=1/h^4.*spdiags(ones(Sys_dim,1)*[1 -4 6 -4 1],-2:2,Sys_dim,Sys_dim);
tmp=full(K_1);
ws(ii,:)=sqrt(eig(K_1));
end
plot(L,ws,[5 15],[1 1],[5 15],[sqrt(3) sqrt(3)])