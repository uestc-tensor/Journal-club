function [A,r]=trbals(T,maxiter,tot)
% appriximate origin data within certain error rho.
%
%--input:
%       T: origin tensor
%       rho: desired approximation errors
%       maxiter: max iteration number
%--output:
%       A: tensor ring decomposition of origin tensor T
%       r: tensor ring ranks
%---------------------------
dim=size(T);
d=size(dim);
d=d(2);
tnorm=norm(reshape(T,1,[]),'fro');
sigma=zeros(1,d);
for k=1:d
    if k == 1
        sigma(1,k)=(sqrt(2)*tot*tnorm)/sqrt(d);
    else
        sigma(1,k)=(tot*tnorm)/sqrt(d);
    end
end

%initialization
r=ones(d+1,1);
A=cell(1,d);
for k=1:d
    %     r(k),r(k+1)=1;
    if k==d
        A{k}=rand(r(k),dim(k),r(1));
    else
        A{k}=rand(r(k),dim(k),r(k+1));
    end
end
% r(d+1)=r(1);
l=0;
%loop
rho=tot;
% Amatrix=cell(1,d);
while l<maxiter
    l=l+1;
    for k=1:d
        
        %compute the subchain and matrix
        [reg]=updatebreg(A,k,r);
        %upgate Z^(k,k+1)
        if k==d
            Tk=reshape(permute(T,[d,1:1:d-1]),dim(d)*dim(1),[]);
            Amatrix=Tk/reg;
            Amatrix=permute(reshape(Amatrix,[dim(d)*dim(1),r(d),r(2)]),[2,1,3]);
            Amatrix=reshape(Amatrix,[r(d)*dim(d),dim(1)*r(2)]);
            %% SVD
            %         sigma=max(tot*tnorm/sqrt(d),rho*tnorm/sqrt(d));
            [U,S,V,r(1)]=SVT(Amatrix, sigma(1,k) );
            A{k}=reshape(U,[r(d),dim(d),r(1)]);
            A{1}=reshape(S*V',[r(1),dim(1),r(2)]);
            r(d+1)=r(1);
        else
            Tk=reshape(permute(T,[k,k+1:1:d,1:1:k-1]),dim(k)*dim(k+1),[]);
            Amatrix=Tk/reg;
            Amatrix=permute(reshape(Amatrix,[dim(k)*dim(k+1),r(k),r(k+2)]),[2,1,3]);
            Amatrix=reshape(Amatrix,[r(k)*dim(k),dim(k+1)*r(k+2)]);
            %% SVD
            %          sigma=max(tot*tnorm/sqrt(d),rho*tnorm/sqrt(d));
            [U,S,V,r(k+1)]=SVT(Amatrix, sigma(1,k) );
            A{k}=reshape(U,[r(k),dim(k),r(k+1)]);
            A{k+1}=reshape(S*V',[r(k+1),dim(k+1),r(k+2)]);
        end
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