function y = expm1div(x)
% computes (exp(x)-1)/x
% when x is a float or an interval 
% and x can be either real or complex

if altisintval(x(1))
    if isreal(x(1))
        % real intervals case
        x1=intval(inf(x));
        x2=intval(sup(x));
        x1(x1>1000)=intval(1000); % taking x1 larger does not improve the evaluation
        x2(sup(x)>realmax)=infsup(realmax,Inf); % catch the edge case x2=Inf
    
        % error estimate 
        % sum_k>=3 r^k/(k+1)! <= r^3/24*sum_k>=0 4!/(k+4)! r^k <= r^3/24*sum_k>=0 k! r^k
        w1=1+x1/2+x1.^2/6-abs(x1).^3/24.*exp(abs(x1));
        w2=1+x2/2+x2.^2/6+abs(x2).^3/24.*exp(abs(x2));
    
        v1=(exp(x1)-1)./x1;
        v1(x1==0)=1;
        v2=(exp(x2)-1)./x2;
        v2(x2==0)=1;
        
        y1=max(inf(w1),inf(v1));
        y2=min(sup(w2),sup(v2));
        y=infsup(y1,y2);
    else
    % complex case x=z as midrad
        zm=intval(mid(x));
        zr=intval(rad(x));
        y1m=1+zm/2+zm.^2/6;
        y1r1=intval(rad(y1m));
        y1m=complex(mid(y1m));
     
        % error comestimate see above    
        y1r2=abs(zm).^3/24.*exp(abs(zm));
        % error estimate sup_z(|f'(z)|)*zr <= sup_z(f(|z|)*zr
        y1r3=sup(expm1div(abs(x))).*zr;
        y1r=sup(y1r1+y1r2+y1r3);

        y1=midrad(y1m,y1r);
        y=y1;
        % when ball does not contain 0
        y2=(exp(x)-1)./x;
        nonnan=~isnan(y2);
        y(nonnan)=intersect(y1(nonnan),y2(nonnan));
    end
else
    % floats    
    y=expm1(x)./x;
    y(x==0)=1;
end

