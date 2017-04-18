function testres=btss(k,n,t,boots_naux,target)
if nargin<5, target=1:length(t); end
boots=boots_naux.boots;
naux=boots_naux.naux;
groups=numel(k);
bootnum=size(boots{1},1);
for i=1:groups
    K(i,:)=wa(k{i},n{i});
    N(i,1)=sum(n{i});
end
K_null=wa(K,N);
w=power(t,-2);
if t(1)==0, w(1)=0; end
for i=1:groups
    btss(i,:)=w.*power(K(i,target)-K_null(target),2);
end
BTSS=sum(N'*btss);clear btss
for i=1:bootnum
    for j=1:groups
        btss(j,:)=w.*power(boots{j}(i,target)-K_null(target),2);
        N(j,1)=sum(naux{i}{j});
    end
    BTSSboot(i)=sum(N'*btss);
end
[p xBTSS fBTSS]=pval(BTSS,BTSSboot,[]);
testres.p=1-fBTSS(length(fBTSS)-(nnz(BTSS<xBTSS)+...
    (nnz(BTSS<xBTSS)==length(fBTSS))*(length(fBTSS)-1)));
testres.BTSS=BTSS;
testres.BTSSdist=[xBTSS;fBTSS];
end