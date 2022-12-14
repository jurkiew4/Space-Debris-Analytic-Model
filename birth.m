function value=birth(r)



c=[200 850 500 1500 700];
k=[.00001, .25, .25, 0, .15];
w=[7.07 10 20 100 100];

coeff=k(1)*w(1)*sqrt(pi)*(w(1)^2+(6378+200)^2)+k(2)*w(2)*sqrt(pi)*(w(2)^2+(6378+850)^2)+k(3)*w(3)*sqrt(pi)*(w(3)^2+(500+6378)^2)+k(4)*w(4)*sqrt(pi)*(w(4)^2+(1500+6378)^2)+k(5)*w(5)*sqrt(pi)*(w(5)^2+(700+6378)^2);
coeff=1/coeff/9.0199;
value=coeff*(k(1)*exp(-(r-200-6378).^2/w(1)^2)+k(2)*exp(-(r-850-6378).^2/w(2)^2)+k(3)*exp(-(r-500-6378).^2/w(3)^2)+k(4)*exp(-(r-1500-6378).^2/w(4)^2)+k(5)*exp(-(r-700-6378).^2/w(5)^2));
end
