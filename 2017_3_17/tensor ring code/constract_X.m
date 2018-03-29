function [al]=constract_X(core,r,d,n)
% constract origin tensor with tensor ring format


    al= core{1};
    A=n(1);
    for i=2:d-1
       cr=core{i};
       cr=reshape(cr,[r(i),n(i)*r(i+1)]);
       
       al=reshape(al,[r(1)*(prod(A)),r(i)]);
       al=al*cr;
       A=[A;n(i)];
    end
    
    cr=core{d};
    cr=permute(cr,[3,1,2]);
    cr=reshape(cr,[r(d)*r(1),n(d)]);
    al=reshape(al,[r(1),numel(al)/(r(1)*r(d)),r(d)]);
    al=permute(al,[2,1,3]);
    al=reshape(al,[numel(al)/(r(1)*r(d)),r(1)*r(d)]); 
    al=al*cr;
end

