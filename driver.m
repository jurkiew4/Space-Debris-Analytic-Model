clf
scale=10;
tbound=365*10*scale;
rmesh=750;
tmesh=10000;
r=linspace(200,2000, rmesh)+6378;
t=linspace(0,tbound, tmesh);
sampst=100*scale*linspace(1,10,10);
sampsr=2*linspace(1,375,375);
U=zeros(1000,rmesh);
ut=zeros(1,rmesh);
ur=ut;
urr=ut;
diff=ut;
coll=ut;
Maxes=zeros(1,1000);
M=5.972e24;
G=6.674e-11*(86400)^2/1000^3;
beta=.1/1000^2;
gammaBar=2000;
Re=6378;
D=diffuse(r);
epoch=2022;

y=[121 86 96 88 74 72 95 111 109 125 120 129 134 210 241 222 221 456 453 586 1274 1809 1980];
tspan=2000:2022;
tchip=pchip(tspan,y/365);




bc1=0;
bc2=0;


hi=1800/length(r);
hj=tbound/tmesh;

u=initial(r);
u(1)=0;
u(end)=u(end-1);
U(1,:)=u;

   upp=pchip(r,u);
   myfun=@(x) ppval(upp,x).*4.*pi.*x.^2;
   Maxes(1)=integral(myfun,r(1),r(end));

for j=2:tmesh
   u1=zeros(1,length(ut));
   u2=u1;

    try 
   ur=gradient(u,hi);
   diff(1)=1/(r(1)^2*hi)*(D(2)*r(2)^2*ur(2)-D(1)*r(1)^2*ur(1));
   for i=2:rmesh-1
       diff(i)=1/(r(1)^2*2*hi)*(D(i+1)*r(i+1)^2*ur(i+1)-D(i-1)*r(i-1)^2*ur(i-1));
   end
   diff(end)=1/(r(end)^2*hi)*(D(end)*r(end)^2*ur(end)-D(end-1)*r(end-1)^2*ur(end-1));
   
   for i=1:rmesh
       coll(i)=1/sqrt(2)*beta*gammaBar*sqrt(G*M/(r(i)))*u(i)^2;
   end
   
   if(epoch+j*hj/365<2022)
    T=ppval(tchip, t);
    elseif(epoch+j*hj/365>=2022)
    T=y(end)/365;
    end
   %Delta=T*birth(r);
   %Delta=-.01/365*u;
   %Delta=birth(r).*exp(-t(i)/100);
   Delta=0;
   ut=diff+coll+Delta;
 
   
   
   u1(2:end-1)=u(2:end-1)+hj*ut(2:end-1);
   u1(1)=0;
   u1(end)=u(end-1)+hi*bc2;
   
   CrankFun=@(x)10000*(model(r,x,hj,epoch+j*hj/365)-u);
   u2=fsolve(CrankFun,u);
   u2(1)=0;
   u2(end)=u2(end-1)+hi*bc2;
   
   u=.5*(u1+u2);
   
   
%    u=.5*(u1+u2);
   
   display(int2str(j))
   
   %plot(r,u)
   %ylim([0,5e-7])
   %pause(.05)
   
%    if(ismember(j,sampst))
%      l=find(sampst==j);
%      upp=pchip(r,u);
%      myfun=@(x) ppval(upp,x).*4.*pi.*x.^2;
%      Maxes(l)=integral(myfun,r(1),r(end));
     U(j,:)=u;
     
%    end
      
   catch 
       plot(hj/365*1000*scale*(1:1000), Maxes(1:1000))
       xlim([0, 10*scale])
       ylim([0,50000])
       display('Something is wrong...')
       break
    end
 
   
end
