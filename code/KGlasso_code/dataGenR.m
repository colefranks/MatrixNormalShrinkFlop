function data = dataGenR(p,f,n)

data = [];

for m=1:n,
    W = randn(f,p);
    data = [data W(:)];
end
clear W;
