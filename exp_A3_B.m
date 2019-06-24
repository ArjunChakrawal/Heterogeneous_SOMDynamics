function e=exp_A3_B(A,B)
    e=1/(numel(A)-1)*sum(sum(( (A-mean2(A)).^3) .* (B-mean2(B))));
    
    