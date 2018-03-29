
function [A,r]=trsvd(T,rho)
% appriximate origin data within certain error rho.
%  
%--input:
%       T: origin tensor
%       rho: desired approximation errors
%--output:
%       A: tensor ring decomposition of origin tensor T
%       r: tensor ring ranks
%---------------------------
dim=size(T);
d=size(dim);
d=d(2);
tnorm=norm(reshape(T,1,[]),'fro');

%compute truncation threshold sigma(k) for k>1.
sigma=zeros(1,d);
for k=1:d
    if k == 1
        sigma(1,k)=(sqrt(2)*rho*tnorm)/sqrt(d);
    else
        sigma(1,k)=(rho*tnorm)/sqrt(d);
    end
end

%initialization
T1=reshape(T,dim(1),prod(dim(2:d)));
[U,S,V,m]=SVT(T1, sigma(1,1) );
%obtain r1,r2
M=round(sqrt(m));
while 1
    b=m/M;
    if b-round(b)==0
        r(1)=M;
        r(2)=b;
        break;
    else
        M=M-1;
    end
end
A =cell(1,d);
B=cell(1,d);
A{1}=permute(reshape(U,[dim(1),r(1),r(2)]),[2,1,3]);
B{1}=permute(reshape(S*V',[r(1),r(2),prod(dim(2:d))]),[2,3,1]);

%loop

for k=2:d-1
    B{k-1}=reshape(B{k-1},[r(k)*dim(k),prod(dim(k+1:d))*r(1)]);
    [U,S,V,r(k+1)]=SVT(B{k-1}, sigma(1,k) );
    A{k}=reshape(U,[r(k),dim(k),r(k+1)]);
    B{k}=reshape(S*V',[r(k+1),prod(dim(k+1:d)),r(1)]);
end
A{d}=B{k};
end
   