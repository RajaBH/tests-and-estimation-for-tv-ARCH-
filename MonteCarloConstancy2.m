function y = MonteCarloConstancy2(STAT,T,b)
%Il s'agit de calculer les p-valeurs pour tester la constance de chaque
%coefficient d'un ARCH(p) (b=largeur de fenêtre, T=taille échantillon)
p=length(STAT)-1;
nbsimu=2000;bopt=b;
STATMC=zeros(nbsimu,p+1);
for s=1:nbsimu
    s
    U=randn(1,T).^2;
    %Calculs des poids P.
chi=zeros(p+1,T);P=zeros(1,T);PchiU=zeros(p+1,T);
Pchichi=zeros(p+1,p+1,T);
chi(1,p+1:T)=1;
for i=2:p+1
    chi(i,p+1:T)=U(p+1-i+1:T-i+1);
end;    
g=sum(chi(:,p+1:T),1);P(p+1:T)=1./(g.^2);   

w=weight(bopt,T);w=w(:,p+1:T);

%Calcul des estimateurs non paramétriques le long de la grille
%ainsi que de la matrice de variance O. 

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

sigma2chap=sum(ahat(:,p+1:T).*chi(:,p+1:T),1);
Inter=(P(p+1:T).^2).*((U(p+1:T)-sigma2chap).^2);

InterPchichi=zeros(p+1,p+1,T);
for i=1:p+1
        for j=1:p+1
        InterPchichi(i,j,:)=w*((Inter.*chi(i,p+1:T).*chi(j,p+1:T))');
        end;
end;  
    centrage=zeros(1,T);renorm=zeros(1,T);O=zeros(T,p+1);
for t=1:T
    M=InterPchichi(:,:,t);
    S=Pchichi(:,:,t);
    O(t,:)=diag((S^(-1))*M*(S^(-1)));
    
end;


%estimateurs sous l'hypothèse nulle.
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
    
    
    
    
    
    
    
for j=1:p+1
STATMC(s,j)=(1/T)*(ahat(j,:)-aparam(j))*((ahat(j,:)-aparam(j))');
varpi1=(1/T)*sum(O(:,j));varpi2=(1/T)*sum(O(:,j).^2);
STATMC(s,j)=T*sqrt(b)*(STATMC(s,j)-(varpi1*(9/15)/(T*b)));
STATMC(s,j)=STATMC(s,j)/(2*sqrt(0.2169*varpi2));
end;   

end;
pval=zeros(1,p+1);
for j=1:p+1
    pval(j)=mean((STATMC(:,j)>STAT(j)));
end;
y=pval;
end

