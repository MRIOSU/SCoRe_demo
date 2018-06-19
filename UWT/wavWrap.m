function [U, Ut, ind] = wavWrap(x, p)

dict = p.dict;
n    = p.n;
N    = p.N; % level of decomposition

%% For 3D imaging undecimated wavelet
nddwt = cell(1,numel(dict));
for i = 1:numel(dict)
    nddwt{i} = nd_dwt_3D(dict{i},n);
end
if numel(dict)==1
    U  = @(x) reshape(nddwt{1}.dec(x,N),[],1);
    Ut = @(x) nddwt{1}.rec(x,N);
end
sz = prod(n);
ind = [0,sz:sz:(7*N+1)*sz*numel(dict)];
