function y = weightCVP(b,T,p)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A=zeros(T,T);
for i=1:T-p
           A(i,:)=i-[1:T];
           A(i,i:i+p)=10*T;
end;
for j=1:p
    A(T-p+j,1:T)=T-p+j-[1:T];A(T-p+j,T-p+j:T)=10*T;

end;
A=A/T;
y=0.75*(1/b)*(1-((A/b).^2)).*(abs(A)<=b);

end

