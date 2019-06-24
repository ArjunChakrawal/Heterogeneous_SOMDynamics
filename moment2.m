function e=moment2(A,n)
    e=1/(numel(A)-1) * sum(sum((A-mean2(A)).^n));
    
    