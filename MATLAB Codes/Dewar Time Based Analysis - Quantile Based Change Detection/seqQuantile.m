function [seqQuantEst,n,m] = seqQuantile(res,seqQuantEst,n,m,tau,N)

    s = 0.02;
  
        if(res >= seqQuantEst)
           n = n + 1;
        else
           m = m + 1;
        end

        if(n > N*tau)
           seqQuantEst = seqQuantEst + s;
           n = 0;
           m = 0;

        elseif(m > N*(1 - tau))
           seqQuantEst = seqQuantEst - s;     
           n = 0;
           m = 0;

        elseif(m + n == N)
           n = 0;
           m = 0;
        end
   
end