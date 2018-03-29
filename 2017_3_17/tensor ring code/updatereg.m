function [a] = updatereg(core, ind,dim,r)
% update subchain Z^(~=ind)
%  
%--input:
%       core: cores of the tensor ring decomposition 
%       ind: the index of current subchain
%       r: tensor ring ranks
%--output:
%       a: subchain Z^(~=ind)
%---------------------------

 d=length(core); 
 n=dim;  

%     T=permute(core{ind},[k,k+1:1:d,1:1:k-1]),dim(k),[]);
%% compute reg for ind=d
if ind==d
    if ind-1==1
        al = reshape(permute(reshape(core{1},[r(1),n(1),r(2)]),[2,3,1]),[n(1),(r(1)*r(2))]) ;
    else
        al= core{1};
        for i=2:ind-1
           cr=core{i};
           cr=reshape(cr,[r(i),n(i)*r(i+1)]);
           al=reshape(al,[numel(al)/r(i),r(i)]);
           al=al*cr;
        end
        al=reshape(al,[r(1),numel(al)/(r(ind)*r(1)),r(ind)]); 
    end
%     a=al;
    a=reshape(permute(al,[3,1,2]),[(r(ind)*r(1)),numel(al)/(r(ind)*r(1))]);
    return;
end

%% comput reg for ind=1
if ind==1
    if ind+1==d
        ar=reshape(core{d},[r(d)*n(d),r(d+1)]);
    else
        ar=core{d};
        for i=d-1:-1:ind+1
          cr=core{i};
          cr=reshape(cr,[r(i)*n(i),r(i+1)]);
          ar=reshape(ar,[r(i+1),numel(ar)/r(i+1)]);
          ar=cr*ar;
        end
        ar=reshape(ar,[r(ind+1),numel(ar)/(r(ind+1)*r(1)),r(1)]);  
    end
%     a=ar;
    a=reshape(permute(ar,[3,1,2]),[(r(ind+1)*r(1)),numel(ar)/(r(ind+1)*r(1))]);
    return;
end
 %% compute reg for ind~=1 and ind~=d  
 %compute left factors
if ind-1==1
    al = reshape(core{1},[r(1),n(1)*r(2)]) ;
else
    al= core{1};
    for i=2:ind-1
       cr=core{i};
       cr=reshape(cr,[r(i),n(i)*r(i+1)]);
       al=reshape(al,[numel(al)/r(i),r(i)]);
       al=al*cr;
    end
%     al=reshape(al,[numel(al)/r(ind-1),r(ind-1)]); 
%  al=reshape(al,[r(1),numel(al)/(r(1)*r(ind-1)),r(ind-1)]);  
 al=reshape(al, r(1),numel(al)/r(1));
end

%compute right factors
if ind+1==d
    ar=reshape(core{d},[n(d)*r(d),r(1)]);
else
    ar=core{d};
    for i=d-1:-1:ind+1
      cr=core{i};
      cr=reshape(cr,[r(i)*n(i),r(i+1)]);
      ar=reshape(ar,[r(i+1),numel(ar)/r(i+1)]);
      ar=cr*ar;
    end
    ar=reshape(ar,[numel(ar)/r(1),r(1)]);  
end
a=ar*al;
a=reshape(a,[r(ind+1),numel(a)/(r(ind+1)*r(ind)),r(ind)]);
a=permute(a,[3,1,2]);
a=reshape(a,[(r(ind+1)*r(ind)),numel(a)/(r(ind+1)*r(ind))]);

% a=reshape(ar,[numel(ar)/r(1),r(1)])*reshape(al,[r(1),numel(al)/r(1)]);
% a=reshape(a,[r(ind+1),numel(a)/(r(ind+1)*r(ind)),r(ind)]);
% a=permute(a,[2,3,1]);
% a=reshape(a,[numel(a)/(r(ind+1)*r(ind)),(r(ind+1)*r(ind))]);

return;
