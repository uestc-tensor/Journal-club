function [a] = updatebreg(core,ind,r)
% update subchain Z^(~=ind,ind+1)
%  
%--input:
%       core: cores of the tensor ring decomposition 
%       ind: the index of current subchain
%       r: tensor ring ranks
%--output:
%       a: subchain Z^(~=ind,ind+1)
%---------------------------

d=length(core); 

l=0;
if ind==d
    index=[2:1:ind-1];
    rleft=r(2);
else
    index=[ind+2:1:d,1:1:ind-1];
    rleft=r(ind+2);
end
Z=cell(1,length(index));
for i=index
    l=l+1;
    Z{l}=core{i};
end


rright=r(ind);
r=r(index);
r(1)=1;
r(length(index)+1)=1;

a= Z{1};
for i=2:length(index)
       cr=Z{i};
       cr=reshape(cr,r(i),[]); 
       a=reshape(a,[],r(i));
       a=a*cr;
end

a=reshape(a,rleft,[],rright);
a=permute(a,[3,1,2]);
a=reshape(a,(rleft*rright),[]);



return;
