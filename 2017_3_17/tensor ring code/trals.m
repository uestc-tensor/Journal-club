function [A]=trals(T,r,maxiter,tot)
% appriximate origin data within certain error rho.
%
%--input:
%       T: origin tensor
%       r: given tensor ring ranks
%       rho: desired approximation errors
%       maxiter: max iteration number
%--output:
%       A: tensor ring decomposition of origin tensor T
%---------------------------
dim=size(T);
d=size(dim);
d=d(2);
tnorm=norm(reshape(T,1,[]),'fro');
%initialization of tensor ring cores
A=cell(1,d);
for k=1:d
    if k==d
        A{k}=randn(r(k),dim(k),r(1));
    else
        A{k}=randn(r(k),dim(k),r(k+1));
    end
end
r(d+1)=r(1);
l=0;
%loop

while l<maxiter
    l=l+1;
    for k=1:d
        
        %compute the subchain and matrix
        [reg]=updatereg(A,k,dim,r);
        %upgate Z
        Tk=reshape(permute(T,[k,k+1:1:d,1:1:k-1]),dim(k),[]);
        Amatrix=Tk/reg;
        %nomalize
        %         if k~=d
        %             Amatrix=normc(Amatrix);
        %         end
        A{k}=permute(reshape(Amatrix,[dim(k),r(k),r(k+1)]),[2,1,3]);
        
        %relative error
        T_pre=constract_X(A,r,d,dim);
        T_pre=reshape(T_pre,dim);
        rho=norm(reshape(T_pre-T,1,[]),'fro');
        rho=rho/tnorm
        
    end
    
    if rho<tot
        break
    end
end