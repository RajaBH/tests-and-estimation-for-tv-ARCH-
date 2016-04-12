function y=TestTVcoefs(X,p)

%X is a row vector. p the number of lags.
%This function gives the p-values for testing 
%the hypotheses: a_j constant for j=1...p
%and the hypothesis: (a_1,...,a_p) constant.

U=X.^2;
T=length(U);
mup=mean(U);

bopt=CVbis(X,p);

chi=zeros(p+1,T);P=zeros(1,T);PchiU=zeros(p+1,T);
Pchichi=zeros(p+1,p+1,T);
chi(1,p+1:T)=mup;
for i=2:p+1
    chi(i,p+1:T)=U(p+1-i+1:T-i+1);
end;    
g=sum(chi(:,p+1:T),1);P(p+1:T)=1./(g.^2);   
chi(1,p+1:T)=1;
w=weight(bopt,T);w=w(:,p+1:T);

ahat=zeros(p+1,T);

for i=1:p+1
        PchiU(i,:)=w*((P(p+1:T).*U(p+1:T).*chi(i,p+1:T))');
        for j=1:p+1
        Pchichi(i,j,:)=w*((P(p+1:T).*chi(i,p+1:T).*chi(j,p+1:T))');
        end;
end;  
for t=1:T   
       ahat(:,t)=Pchichi(:,:,t)\PchiU(:,t);
end;
ahattronc=max(ahat,0);
sigma2chap=sum(ahattronc(:,p+1:T).*chi(:,p+1:T),1);

Inter=(P(p+1:T).^2).*((U(p+1:T)-sigma2chap).^2);

InterPchichi=zeros(p+1,p+1,T);
for i=1:p+1
        for j=1:p+1
        InterPchichi(i,j,:)=w*((Inter.*chi(i,p+1:T).*chi(j,p+1:T))');
        end;
end;  
    centrage=zeros(1,T);renorm=zeros(1,T);O=zeros(T,p+1);
    centrage2=zeros(1,T);renorm2=zeros(1,T);
for t=1:T
    M=InterPchichi(:,:,t);
    S=Pchichi(:,:,t);
    mamat=(S^(-1))*M*(S^(-1));
    O(t,:)=diag(mamat);
    centrage2(t)=trace(mamat);
    renorm2(t)=trace(mamat^2);
    centrage(t)=trace(mamat(2:p+1,2:p+1));
    renorm(t)=trace(mamat(2:p+1,2:p+1)^2);
end;


aparam=zeros(1,p+1);Q1=zeros(p,T,p+1);Q2=zeros(p,T,p+1);

for j=1:p+1
    PJJ=Pchichi([1:j-1,j+1:p+1],[1:j-1,j+1:p+1],:);
    PJU=PchiU([1:j-1,j+1:p+1],:); 
    PJK=Pchichi([1:j-1,j+1:p+1],j,:);PJK=permute(PJK,[1,3,2]);

   for t=p+1:T
    Q1(:,t,j)=PJJ(:,:,t)\PJU(:,t);
    Q2(:,t,j)=PJJ(:,:,t)\PJK(:,t);
   end;
end;

for j=1:p+1
    V=U(p+1:T)-sum(chi([1:j-1,j+1:p+1],p+1:T).*Q1(:,p+1:T,j),1);
    W=chi(j,p+1:T)-sum(Q2(:,p+1:T,j).*chi([1:j-1,j+1:p+1],p+1:T),1);
    aparam(j)=(sum(P(p+1:T).*V.*W))/sum(P(p+1:T).*(W.^2));
end;

   S3=P(p+1:T)*(w');
    Q1=(P(p+1:T).*U(p+1:T))*(w');
    Q1=Q1./S3;

    Q2=zeros(p,1,T);K=zeros(p,1,T);
    S3=repmat(S3,[p,1,p]);
    S3=permute(S3,[1,3,2]);
for i=1:p
    Q2(i,1,:)=w*((P(p+1:T).*U(p+1-i:T-i))');
    K(i,1,p+1:T)=U(p-i+1:T-i);
end;
Q2=repmat(Q2,[1,p,1]);Q2=Q2./S3;
K=repmat(K,[1,p,1]);
diff=K(:,:,p+1:T)-Q2(:,:,p+1:T);difft=permute(diff,[2,1,3]);
fact=repmat(P(p+1:T),[p,1,p]);
fact=permute(fact,[1,3,2]);
D=mean(fact.*diff.*difft,3);
diff2=P(p+1:T).*(U(p+1:T)-Q1(p+1:T));
diff2=repmat(diff2,[p,1,1]);
diff2=permute(diff2,[1,3,2]);
N=mean(diff2.*diff(:,1,:),3);
N=permute(N,[1,3,2]);
aparam2=D\N;


STAT=zeros(1,p+2);
for j=1:p+1
STAT(j)=(1/T)*(ahat(j,:)-aparam(j))*((ahat(j,:)-aparam(j))');
varpi1=(1/T)*sum(O(:,j));varpi2=(1/T)*sum(O(:,j).^2);
STAT(j)=T*sqrt(bopt)*(STAT(j)-(varpi1*(9/15)/(T*bopt)));
STAT(j)=STAT(j)/(2*sqrt(0.2169*varpi2));
end;
difference=ahat(2:p+1,:)-repmat(aparam2,[1,T]);
STAT(p+2)=(1/T)*sum(sum(difference.^2));
varpi1=mean(centrage);
varpi2=mean(renorm);
STAT(p+2)=T*sqrt(bopt)*(STAT(p+2)-(varpi1*(9/15)/(T*bopt)));
STAT(p+2)=STAT(p+2)/(2*sqrt(0.2169*varpi2));


pvaleur=MonteCarloConstancy2(STAT,T,bopt);
y=pvaleur;
