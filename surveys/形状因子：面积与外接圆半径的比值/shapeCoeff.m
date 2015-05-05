function [coeff,R,S] = shapeCoeff( a,b,c )
% Shape coeff is S/(pi*R^2)
    p=(a+b+c)/2;
    S=sqrt(p.*(p-a).*(p-b).*(p-c));
    R=a.*b.*c/4./S;
    coeff=S./(pi*R.^2);
end

