function value=diffuse(r)
alpha=.578329;
dend=.0001;

lambda=log(alpha/dend)/1000;

% lambda=.0134;
% alpha=69.168;
mid=zeros(1,length(r));
for i=1:length(r)
if(r(i)<1000+6378)
    mid(i)=alpha*exp(-lambda*(r(i)-6378));
elseif(r(i)>=1000+6378)
    mid(i)=dend;
end
end
value=mid;
end
