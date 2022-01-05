function x_d =eq_of_var(t,x,Om,phi,pars)

% k=pars(1);
% c=pars(2);
% ampl(1)=pars(3);
% ampl(2)=pars(4);
% -pars(4)*cos(Om(2)*t+phi(2))

A=[0 1; -pars(1)-pars(3)*cos(Om(1)*t+phi(1))-pars(4)*cos(Om(2)*t+phi(2)) -pars(2)];

tmp=A*([x(1:2)  x(3:4) ]);%([x(1:2).';  x(3:4).'  ]);
x_d=tmp(:);%reshape(tmp.',4,1);%

    
end