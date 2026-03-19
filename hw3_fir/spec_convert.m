function [A,B] = spec_convert(C,D,typein,typeout)

% Step 1: convert input → rel (Ap, As)
switch typein
    
    case 'rel'
        Ap = C;
        As = D;
        
    case 'abs'
        dp = C;
        ds = D;
        
        Ap = 20*log10((1+dp)/(1-dp));
        As = 20*log10((1+dp)/ds);
        
    case 'ana'
        eps = C;
        A_ = D;
        
        Ap = 10*log10(1 + eps^2);
        As = 20*log10(A_);
end


% Step 2: convert rel → output
switch typeout
    
    case 'rel'
        A = Ap;
        B = As;
        
    case 'abs'
        dp = (10^(Ap/20)-1)/(10^(Ap/20)+1);
        ds = (1+dp)/10^(As/20);
        
        A = dp;
        B = ds;
        
    case 'ana'
        eps = sqrt(10^(0.1*Ap) - 1);
        A_ = 10^(0.05*As);
        
        A = eps;
        B = A_;
end

end