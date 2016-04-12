function y=EstimSParch22(X,I,p)


U=X.^2;
T=length(U);
sup=length(I);
K=zeros(p,1,T);


pond=ones(1,T);
for i=1:p
    K(i,1,p+1:T)=U(p-i+1:T-i);
end;
 K=repmat(K,[1,p,1]);

K2=K(:,1,:);K2=permute(K2,[1,3,2]);
%Cross-validation.
prev=zeros(1,sup);
for u=1:sup
    b=I(u); 
    w=weightCVP(b,T,p);w=w(:,p+1:T);  
    
    
a0chap=zeros(1,T);Q2=zeros(p,1,T);

S3=sum(w',1);
Q1=U(p+1:T)*(w');
Q1=Q1./S3;

S3=repmat(S3,[p,1,p]);
S3=permute(S3,[1,3,2]);
for i=1:p
    Q2(i,1,:)=w*(U(p+1-i:T-i)');
end;
Q2=repmat(Q2,[1,p,1]);Q2=Q2./S3;
diff=K(:,:,p+1:T)-Q2(:,:,p+1:T);difft=permute(diff,[2,1,3]);
fact=repmat(pond(p+1:T),[p,1,p]);
fact=permute(fact,[1,3,2]);
D=mean(fact.*diff.*difft,3);
diff2=pond(p+1:T).*(U(p+1:T)-Q1(p+1:T));
diff2=repmat(diff2,[p,1,1]);
diff2=permute(diff2,[1,3,2]);
N=mean(diff2.*diff(:,1,:),3);
N=permute(N,[1,3,2]);
aparam=D\N;
Q2=permute(Q2(:,1,:),[1,3,2]);
a0chap(p+1:T)=Q1(p+1:T)-(aparam')*Q2(:,p+1:T);
 prev(u)=mean(pond(p+1:T).*((U(p+1:T)-a0chap(p+1:T)-(aparam')*K2(:,p+1:T)).^2));   
end;

[opt,indopt]=min(prev);bopt=I(indopt);




y=bopt;