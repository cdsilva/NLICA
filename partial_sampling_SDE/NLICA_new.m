function [V_NLICA, D_NLICA] = NLICA_new(data, inv_c, ep, N)

% [m, n] = size(data);
% 
% W = zeros(m, m);
% 
% for i=1:m
%     W(:,i) = sum((ones(m,1) * data(i,:) - data) * squeeze(inv_c(i,:,:)) ...
%              .* (ones(m,1) * data(i,:) - data),2);
% 
%     %waitbar(i/m)
%     %for j=1:i
%     %    W(i, j) = (data(i,:) - data(j,:)) * squeeze(inv_c(i,:,:) + ...
%     %                                                inv_c(j,:,:)) * ...
%     %              (data(i,:)-data(j,:))';
%     %    W(j,i) = W(i,j);
%     %end
% end
% W = W + W';
% 
% W = exp(-W/(4*eps));
% 
% for i=1:m
%     W(i,:) = W(i,:) / sum(W(i,:));
% end
% 
% [V_NLICA, D_NLICA] = eigs(W, N);
% 
% [~, I] = sort(abs(diag(D_NLICA)), 'descend');
% 
% D_NLICA = D_NLICA(I,I);
% V_NLICA = V_NLICA(:,I);

%%
%inv_c = inv_c_noisy;

%
% Choose subset of examples
%
M = size(data, 1);
m = min(M, 5000);
subidx = 1:m;
dataref = data(subidx, :);

Dis = zeros(M, m);
h = waitbar(0, 'Please wait');
for j = 1:m
    waitbar(j/m, h);
    tmp1 = inv_c(:,:,subidx(j)) * dataref(j,:)';

    a2 = dataref(j,:) * tmp1;
    b2 = sum(data .* (inv_c(:,:,subidx(j)) * data')',2);
    ab = data * tmp1;
    Dis(:,j) = repmat(a2, M, 1) + b2 - 2*ab;
end
close(h);

%%

A = exp(-Dis/(4*ep));
W_sml=A'*A;    
d1=sum(W_sml,1);
A1=A./repmat(sqrt(d1),M,1);
W1=A1'*A1;

d2=sum(W1,1);
A2=A1./repmat(sqrt(d2),M,1);
W2=A2'*A2;

D=diag(sqrt(1./d2));

% computing eigenvectors:
[V,E] = eigs(W2,10);
diag(E(2:end,2:end))
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

V_NLICA = psi_mat;
D_NLICA = E;
