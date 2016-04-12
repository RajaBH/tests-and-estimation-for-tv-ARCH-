function y=CVbis(X,p)

%Cross-validation method for selecting the bandwidth in tv-ARCH 
%processes.


U=X.^2;

    mup=mean(U);
T=length(U);
mini=0.01;maxi=0.3;
pas=0.001;
I=[mini:pas:maxi];
sup=length(I);

prev=zeros(1,sup);
ahat=zeros(p+1,T,sup);
    

chi=zeros(p+1,T);P=zeros(1,T);PchiU=zeros(p+1,T);
Pchichi=zeros(p+1,p+1,T);

chi(1,p+1:T)=mup;
for i=2:p+1
    chi(i,p+1:T)=U(p+1-i+1:T-i+1);
end;    
g=sum(chi(:,p+1:T),1);P(p+1:T)=1./(g.^2);   
chi(1,p+1:T)=1;

for u=1:sup
    w=weightCVP(I(u),T,p);w=w(:,p+1:T);
   
    for i=1:p+1
        PchiU(i,:)=w*((P(p+1:T).*U(p+1:T).*chi(i,p+1:T))');
        for j=1:p+1
        Pchichi(i,j,:)=w*((P(p+1:T).*chi(i,p+1:T).*chi(j,p+1:T))');
        end;
    end;  
    for t=1:T     
       ahat(:,t,u)=Pchichi(:,:,t)\PchiU(:,t);
    
    end;
    past=sum(ahat(:,p+1:T,u).*chi(:,p+1:T),1);
    prev(u)=mean(P(p+1:T).*((U(p+1:T)-past).^2));
end;
[opt,indopt]=min(prev);bopt=I(indopt);
y=bopt;
