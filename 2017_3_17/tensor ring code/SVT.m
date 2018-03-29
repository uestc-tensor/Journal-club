function [U,S,V,r] = SVT(Z, lambda )
% Singular Value Thresholding
% Check if upper bound is lower than threshold first
% if not then proceed with SVT

if ( size(Z,1) > size(Z,2) )
    
    ZZ = Z' * Z;
    
else
    ZZ = Z * Z';
end

if (max( sum( abs( ZZ ), 1 ) ) < lambda^2)
    
    Z = Z * 0;
    
else
    
    [U,S,V] = svd( Z, 'econ' );
    
    Z = U * diag(SoftThresh( diag(S), lambda )) * V';
    S=diag(SoftThresh( diag(S), lambda )) ;
    r=rank(S);
    S=S(1:r,1:r);
    U=U(:,1:r);
    V=V(:,1:r);
end
