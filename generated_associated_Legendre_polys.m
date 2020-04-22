

%% Compute associated Legendre polynomials
%
% $$P_n^m(x) = (-1)^m \cdot 2^n \cdot (1-x^2)^{m/2}\cdot \sum_{k=m}^n \frac{k!}{(k-m)!}\cdot x^{k-m}{n \choose k}{(n+k-1)/2 \choose n}$$
%
% en.wikipedia.org/wiki/Associated_Legendre_polynomials
%%

N = 13;

x = linspace(-1, 1);
Pnm = {};

for n = 1:N
    for m = 0:n
        
    end
end

P = legendre(n, X);
