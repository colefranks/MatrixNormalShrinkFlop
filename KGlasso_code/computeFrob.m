function [out] = computeFrob(A,B,C,D)
%COMPUTEFROB Function computes squared Frobenius error.

p=size(A,1);
out=0;
for i=1:p,
    for j=1:p,
        out = out + norm(A(i,j)*B-C(i,j)*D,'fro')^2;
    end
end


end

