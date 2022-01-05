function x_d =Beam_eq_of_var(t,x,Oms,ampls,C,K_1,K_2,phi0)

 n=min(size(C));
A=[zeros(size(C)) eye(size(C)); -K_1+K_2.*(ampls(1)*cos(Oms(1)*t+phi0(1))+ampls(2)*cos(Oms(2).*t+phi0(2))) -C  ];

%reshape(x,2*n,2*n)-[x(1:2*n)  x(2*n+1:4*n)  x(4*n+1:6*n)  x(6*n+1:8*n)  x(8*n+1:10*n) x(10*n+1:12*n) x(12*n+1:14*n) x(14*n+1:16*n) x(16*n+1:18*n) x(18*n+1:20*n)]
tmp=A*reshape(x,2*n,2*n);%([x(1:n)  x(n+1:2*n) ])
x_d=tmp(:);%reshape(tmp.',4,1);%

    
end