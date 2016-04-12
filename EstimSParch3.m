function y=EstimSParch3(X,p)
%This function computes the estimates of lag coefficients and their s.e
%when lag coefficients are constant.
b=EstimSParch2(X,p);

U=X.^2;
T=length(U);

P=ones(1,T);K=zeros(p,1,T);
for i=1:p
    K(i,1,p+1:T)=U(p-i+1:T-i);
end;
 K=repmat(K,[1,p,1]);
g=sum(K(:,1,p+1:T),1);

mu=mean(U);

P(p+1:T)=1./((mu+g(:)).^2);

w=weight(b,T);w=w(:,p+1:T);
Q2=zeros(p,1,T);

S3=P(p+1:T)*(w');
Q1=(P(p+1:T).*U(p+1:T))*(w');
Q1=Q1./S3;

S3=repmat(S3,[p,1,p]);
S3=permute(S3,[1,3,2]);
for i=1:p
    Q2(i,1,:)=w*((P(p+1:T).*U(p+1-i:T-i))');
end;
Q2=repmat(Q2,[1,p,1]);Q2=Q2./S3;
diff=K(:,:,p+1:T)-Q2(:,:,p+1:T);difft=permute(diff,[2,1,3]);
fact=repmat(P(p+1:T),[p,1,p]);
fact=permute(fact,[1,3,2]);
D=mean(fact.*diff.*difft,3);
diff2=P(p+1:T).*(U(p+1:T)-Q1(p+1:T));
diff2=repmat(diff2,[p,1,1]);
diff2=permute(diff2,[1,3,2]);
N=mean(diff2.*diff(:,1,:),3);
N=permute(N,[1,3,2]);
aparam=D\N;
aparam=max(aparam,0);
Q2=permute(Q2(:,1,:),[1,3,2]);
a0chap=Q1-(aparam')*Q2;
sigma2chap(p+1:T)=a0chap(p+1:T)+(aparam')*(permute(K(:,1,p+1:T),[1,3,2]));

fact1=repmat(P(p+1:T),[p,1,p]);
fact1=permute(fact1,[1,3,2]);
Res=(U(p+1:T)-sigma2chap(p+1:T)).^2;
fact2=repmat((P(p+1:T).^2).*Res,[p,1,p]);
fact2=permute(fact2,[1,3,2]);

Sigma1=mean(fact1.*diff.*difft,3);
Sigma2=mean(fact2.*diff.*difft,3);
Sigma=inv(Sigma1)*Sigma2*inv(Sigma1);
vari=sqrt(diag(Sigma));

y=[aparam',vari'/sqrt(T)];

%y=a0chap;
