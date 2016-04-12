function y=ordertvARCH(X,q)
% Selection of the number of lags. X is a row vector 
% which coordinates are the returns and q is the maximal number of lags
% (q=10 is used in the examples).
U=X.^2;Crit=zeros(1,q+1);
T=length(U);

bopt=CVbis(X,q);     
w=weight(bopt,T);
w=w(:,q+1:T);

Pmax=zeros(1,T);
chimax=zeros(q+1,T);
chimax(1,q+1:T)=mean(U);

for i=2:q+1
    chimax(i,q+1:T)=U(q+1-i+1:T-i+1);
end;    
g=sum(chimax(:,q+1:T),1);Pmax(q+1:T)=1./(g.^2);   
Inter=(Pmax(q+1:T).*U(q+1:T))*(w');
Inter=Inter./(Pmax(q+1:T)*w');
Crit(1)=log(mean(Pmax(q+1:T).*(U(q+1:T)-Inter(q+1:T)).^2))+log(log(T))/(bopt*T);
 
for p=1:q   

chi=zeros(p+1,T);PchiU=zeros(p+1,T);
Pchichi=zeros(p+1,p+1,T);
chi(1,p+1:T)=1;
for i=2:p+1
    chi(i,p+1:T)=U(p+1-i+1:T-i+1);
end;    



ahat=zeros(p+1,T); 
for i=1:p+1
        PchiU(i,:)=w*((Pmax(q+1:T).*U(q+1:T).*chi(i,q+1:T))');
        for j=1:p+1
        Pchichi(i,j,:)=w*((Pmax(q+1:T).*chi(i,q+1:T).*chi(j,q+1:T))');
        end;
end;  
for t=1:T   
       ahat(:,t)=Pchichi(:,:,t)\PchiU(:,t);
end;
ahattronc=max(ahat,0);
sigma2chap=sum(ahattronc(:,q+1:T).*chi(:,q+1:T),1);

Inter=Pmax(q+1:T).*((U(q+1:T)-sigma2chap).^2);
Crit(p+1)=log(mean(Inter))+log(log(T))*(p+1)/(bopt*T);
 end;

[l,n]=min(Crit);
y=n-1;
