function [psi_mat, E] = NL_ICA2(data, inv_c, ep)

[m, n] = size(data);

Dis = zeros(m);
h = waitbar(0, 'Please wait');
for i = 1:m
    waitbar(i/m, h);
    for j = 1:m
        Dis(i,j) = (data(i,:) - data(j,:)) * inv_c(:,:,j) * (data(i,:) - data(j,:))';
    end
end
close(h);

%ep = median(median(Dis));

%%
A = exp(-Dis/(4*ep));
W_sml=A'*A;    
d1=sum(W_sml,1);
A1=A./repmat(sqrt(d1),m,1);
W1=A1'*A1;

d2=sum(W1,1);
A2=A1./repmat(sqrt(d2),m,1);
W2=A2'*A2;

D=diag(sqrt(1./d2));

% computing eigenvectors:
[V,E] = eigs(W2,10);
E(2:end,2:end)
[srtdE,IE] = sort(sum(E),'descend');
V_srt_clds = D*V(:,IE(1,2:10));

% NOW THE EXTENSION:
psi_mat=[];

omega=sum(A2,2);
A2_nrm=A2./repmat(omega,1,m);

for i=1:size(V_srt_clds,2)
    psi=A2_nrm*V_srt_clds(:,i)./sqrt((srtdE(i+1)));
    psi_mat=[psi_mat,psi];
end
