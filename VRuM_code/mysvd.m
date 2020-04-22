function [U,S,V] = mysvd(A)
[n,m] = size(A);
U = eye(n);
V = eye(m);
tol = max(abs(A(:)))*1.e-15;
Arem = inf;

while Arem > tol
  [Qu,Ru] = qr(A);
  U = U*Qu;

    [Qv,Rv] = qr(Ru');
    V = V*Qv;

    A = Rv';
    Arem = norm(tril(A,-1),inf);
  end

U = U(:,1:m);
S = triu(A(1:m,:));

% correct any negative singular values
V = V.*sign(diag(S)).';
S = abs(S);