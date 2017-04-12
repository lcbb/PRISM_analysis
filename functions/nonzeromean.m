function y = nonzeromean(x, dim)
x(x==0) = NaN ;
y = nanmean(x, dim) ;

end