function m = MI(y1, y2, nbins)

A = hist3([y1 y2], [nbins nbins]);

p = A / sum(A(:));

pi = sum(p, 2);
pj = sum(p, 1);


m = 0;
for i=1:nbins
    for j=1:nbins
        if p(i,j) > 0
            m = m + p(i,j)*log2(p(i,j)/(pi(i)*pj(j)));
        end
    end
end

H = 0;
for i=1:nbins
    for j=1:nbins
        if p(i,j) > 0
            H = H - p(i,j) * log2(p(i,j));
        end
    end
end


m = m / H;

