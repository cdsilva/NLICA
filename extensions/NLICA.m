
%subset = [1:2:6];
subset = [1:3];
%subset = [1:6];

inv_c = inv_c_denoised(subset, subset, :);

%data = mean(data, 1);
%data = reshape(data, [size(data, 2), size(data, 3)]);
%data = data';
data = orig_data(:,subset);

%%
%
% Choose subset of examples
%
M = size(data, 1);
m = M;
subidx = 1:m;
dataref = data(subidx, :);

% 
% Anisotropic Kernel
% 
Dis = zeros(M, m);
h = waitbar(0, 'Please wait');
for i = 1:M
    waitbar(i/M, h);
    for j = 1:m
        %Dis(i,j) = [data(i,:) - dataref(j,:)] * [data(i,:) - dataref(j,:)]';
        Dis(i,j) = [data(i,:) - dataref(j,:)] * inv_c(:,:,subidx(j)) * [data(i,:) - dataref(j,:)]';
    end
end
close(h);

%%
%ep = 1e6; % NLICA
ep=1e7;
%ep = 1e6; % DM
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
E(2:end,2:end)
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

%%
figure; 
subplot(1,3,1);
scatter(psi_mat(:,1),psi_mat(:,2),200,orig_data(:,1),'.');
subplot(1,3,2);
scatter(psi_mat(:,1),psi_mat(:,2),200,orig_data(:,2),'.');
subplot(1,3,3);
scatter(psi_mat(:,1),psi_mat(:,2),200,orig_data(:,3),'.');

%%

% Since close to a degenerate case - try to rotate according to [Singer & Coifman, Spectral ICA, 2007]
theta_psi = atan(0.5 * (psi_mat(:, 3).' * (psi_mat(:,1).*psi_mat(:,1) - psi_mat(:,2).*psi_mat(:,2)) ) / ( psi_mat(:,3).' * (psi_mat(:,1) .* psi_mat(:,2)))) / 2
Psi_1 = psi_mat(:, 1) * cos(theta_psi) - psi_mat(:,2) * sin(theta_psi);
Psi_2 = psi_mat(:, 1) * sin(theta_psi) + psi_mat(:,2) * cos(theta_psi);
