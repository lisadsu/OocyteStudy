function out=kboot(k,n,resid,bootnum,type)
%  out=kboot(k,n,resid,bootnum,type)
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 


% default type is within
if nargin<5
    type='within'; 
end

switch type
    case 'within' % bootstrap only one group       
        K_average = wa(k,n);
        N_total = sum(n);
        h_kbootwait = waitbar(0,'Bootstrapping...');
        for i=1:bootnum
            for k=1:length(n)
                if (i/100 - floor(i/100)) == 0
                    try
                        waitbar(i/bootnum, h_kbootwait);
                    catch
                        out = 0;
                        return;
                    end
                end
                b=ceil(length(n)*rand);
                naux(k)=n(b);
                kaux(k,:)=K_average+(1/sqrt(naux(k)))*resid(b,:);
            end
            out(i,:)=naux*kaux(:,:)/sum(naux);
        end
        close(h_kbootwait);
        clear naux kaux
        
    case 'between' % compare between groups
        groups=numel(k);
        for i=1:groups
            K(i,:)=wa(k{i},n{i});
            N(i,1)=sum(n{i});
        end
        for i=1:bootnum
            for j=1:groups
                for k=1:length(n{j})
                    a=ceil(groups*rand);
                    b=ceil(length(n{a})*rand);
                    naux{j}(k)=n{a}(b);
                    kaux{j}(k,:)=wa(K,N)+(1/sqrt(naux{j}(k)))*resid{a}(b,:);
                end
                out.boots{j}(i,:)=naux{j}*kaux{j}(:,:)/sum(naux{j});
                out.naux{i}=naux;
            end
            clear naux kaux
        end
end