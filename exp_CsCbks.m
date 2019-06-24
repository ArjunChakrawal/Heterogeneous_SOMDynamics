function e=exp_CsCbks(A,B,C)
    A=A(:);B=B(:);C=C(:);
    e=(1/(numel(A)-1))*sum((A-mean(A)).*(B-mean(B)).*(C-mean(C)));
    
    