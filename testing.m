for i=1:len_p
    for j=1:len_p
       c(i,j) = myfunc(p(i),y(j)); 
    end
end
contour(c)