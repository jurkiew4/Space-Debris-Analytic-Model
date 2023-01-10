function value=model(r,u, hj,T)
%Initialization of meshes for finite difference approximation
tbound=3;
rmesh=750;
t=linspace(0,tbound, 2000);
ut=zeros(1,rmesh);
hi=1800/750;



%Parametrization of model
M=5.972e24;
G=6.674e-11*(86400)^2/1000^3;
beta=.1/1000^2;
gammaBar=2000;
Re=6378;
D=diffuse(r);


   %Delta=T*birth(r);
   Delta=0;
   ur=gradient(u,hi);
   diff(1)=1/(r(1)^2*hi)*(D(2)*r(2)^2*ur(2)-D(1)*r(1)^2*ur(1));
   for i=2:rmesh-1
       diff(i)=1/(r(1)^2*2*hi)*(D(i+1)*r(i+1)^2*ur(i+1)-D(i-1)*r(i-1)^2*ur(i-1));
   end
   diff(rmesh)=1/(r(end)^2*hi)*(D(end)*r(end)^2*ur(end)-D(end-1)*r(end-1)^2*ur(end-1));
   
   for i=1:rmesh
       coll(i)=1/sqrt(2)*beta*gammaBar*sqrt(G*M/(r(i)))*u(i)^2;
   end
   
value=-hj*(coll+diff+Delta)+u;
end    
