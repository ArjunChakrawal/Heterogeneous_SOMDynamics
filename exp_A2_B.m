function e=exp_A2_B(A,B)
    e=1/(numel(A)-1)*sum(sum(((A-mean2(A)).^2).*(B-mean2(B))));
    
    
 
    
    