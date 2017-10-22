function mylu(n)
A=randn(n,n); b=randn(n,1); Abk=A; pvt = 1:n;
%Factorize A. Your task: transform this part to mydgetrf().
for i = 1 : n-1,
 % pivoting %
 maxind=i; max=abs(A(i,i));
 for t=i+1:n,
 if ( abs(A(t,i))>max )
 maxind=t; max=abs(A(t,i));
 end
 end
 if (max==0)
 disp ( 'LUfactoration failed: coefficient matrix is singular' ); return;
 else
 if (maxind ~= i )
 %save pivoting infomation
 temps=pvt(i);
 pvt(i)=pvt(maxind);
 pvt(maxind)=temps;
 %swap rows
 tempv=A(i,:);
 A(i,:)=A(maxind,:);
 A(maxind,:)=tempv;
 end
 end
 %factorization
 for j = i+1 : n,
 A(j,i) = A(j,i)/A(i,i);
 for k = i+1 : n,
 A(j,k) = A(j,k) - A(j,i) * A(i,k);
 end
 end
end
%verify my factorization with Matlab for small matrix by printing results on screen. May skip in your code
myfactorization=A
mypivoting=pvt
[Matlab_L, Matlab_U, Matlab_P]=lu(Abk)
%forward substitution. Your task: transform this part to mydtrsm().
y(1) = b(pvt(1));
for i = 2 : n,
 y(i) = b(pvt(i)) - sum ( y(1:i-1) .* A(i, 1:i-1) );
end
% back substitution. Your task: transform this part to mydtrsm().
x(n) = y(n) / A(n, n);
for i = n-1 : -1 : 1,
 x(i) = ( y(i) - sum ( x(i+1:n) .* A(i, i+1:n) ) ) / A(i,i);
end
%Matlab solve. Your task: call dgetrf() to factorize and dtrsm() twice (back and forward substit.) to solve.
xx= Abk\b;
%verify my solution with matlab. Your task: verify your solution with the solution from LAPACK.
Solution_Difference_from_Matlab=norm(x'-xx)
