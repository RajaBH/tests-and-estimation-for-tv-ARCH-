function y=TestArcheffectSP(X,p)
%This function gives the p-value for testing the hypothesis
%a1=a2=...=ap=0 in the ARCH model with constant lag coefficients. 
X=Rend;
T=length(X);

bandw=[0.005:0.001:0.05];
bopt=EstimSParch22(X,bandw,p);

valeur=EstimSParch32(X,bopt,p);

Stat=zeros(1,2000);
for j=1:2000
    j 
    Y=randn(1,T);

    Stat(j)=EstimSParch32(Y,bopt,p);
end;

pval=mean((Stat>=valeur));
 y=pval;   
