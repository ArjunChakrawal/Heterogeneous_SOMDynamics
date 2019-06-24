function e=exp_A_B2(A,B)
    e=1/(numel(A)-1)*sum(sum( (A-mean2(A)) .*(B-mean2(B)).^2));
    
    