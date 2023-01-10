function value=initial(r)
A=csvread('DensityNow.csv');
x=linspace(200,2000,750)+6378;
pp=pchip(x,A);


value=ppval(pp, r);

end
