

%
% Choose subset of examples
%
M = size(data, 1);
m = 5000;
subidx = 1:m;
dataref = data(subidx, :);

% 
% Anisotropic Kernel
% 
Dis = zeros(M, m);
% h = waitbar(0, 'Please wait');
% for i = 1:M
%     waitbar(i/M, h);
%     for j = 1:m
%          Dis(i,j) = [data(i,:) - dataref(j,:)] * [data(i,:) - dataref(j,:)]';
%     end
% end
% close(h);

Dis = zeros(M, m);
h = waitbar(0, 'Please wait');
for j = 1:m
    waitbar(j/m, h);

    a2 = dataref(j,:) * dataref(j,:)';
    b2 = sum(data .* data,2);
    ab = data * dataref(j,:)';
    Dis(:,j) = repmat(a2, M, 1) + b2 - 2*ab;
end
close(h);

%ep = 1e1;
%ep = 1e2;
ep = median(median(Dis))

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

Phi_1 = V_srt_clds(:, 1);
Phi_2 = V_srt_clds(:, 2);

% NOW THE EXTENSION:
psi_mat=[];

omega=sum(A2,2);
A2_nrm=A2./repmat(omega,1,m);

for i=1:size(V_srt_clds,2)
    psi=A2_nrm*V_srt_clds(:,i)./sqrt((srtdE(i+1)));
    psi_mat=[psi_mat,psi];
end

%Psi_1 = psi_mat(:,1);
%Psi_2 = psi_mat(:,2);
theta_psi = atan(0.5 * (psi_mat(:, 3).' * (psi_mat(:,1).*psi_mat(:,1) - psi_mat(:,2).*psi_mat(:,2)) ) / ( psi_mat(:,3).' * (psi_mat(:,1) .* psi_mat(:,2)))) / 2
Psi_1 = psi_mat(:, 1) * cos(theta_psi) - psi_mat(:,2) * sin(theta_psi);
Psi_2 = psi_mat(:, 1) * sin(theta_psi) + psi_mat(:,2) * cos(theta_psi);
