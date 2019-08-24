function GammaFit
y = [150  149  151  157  163  169  174  176  177  177  177  181 ...
     185  186  185  184  181  176  173  171  169  167  168  168];
x = 1:24;
p = nlinfit(x,y,@f,[1 10 150])
plot(x,y,'bo',x,f(p,x),'r-')
end

function y = f(abc,x)
a = abc(1); b = abc(2); c = abc(3);
y = c * x.^(a-1) .* exp(-x/b) / (b^a * gamma(a));
end