function y = weight(b,T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A=zeros(T,T);
for i=1:T
           A(i,:)=i-[1:T];
end;

A=A/(T*b);
y=0.75*(1/(T*b))*(1-(A.^2)).*(abs(A)<=1);

end

