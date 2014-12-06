clear all
close all

rng(321);

nsamples = 20;

D = 1;
lambda = linspace(1, 4000, nsamples);
s = sqrt(lambda / D);

neps = 100;
eps_range = logspace(2, 10, neps);

vol_data = zeros(nsamples, neps);

for i=1:nsamples
    W2 = simulate_samples(lambda(i), s(i));
    for j=1:neps
        vol_data(i,j) = sum(exp(-W2(:)/(2*eps_range(j))));
    end
end

%%
figure; 
loglog(eps_range, vol_data)
hold on
loglog(eps_range, 0.01*eps_range)
loglog(eps_range, 0.01*eps_range.^(1/2))

%%

slopes = zeros(nsamples);

% figure;
for i=1:nsamples
    vol1 = vol_data(i,:);
    idx = find(vol1 < (max(vol1)^0.6)*(min(vol1)^0.4) & vol1 > (min(vol1)^0.6)*(max(vol1)^0.4));
%     loglog(eps_range, vol1, '.')
%     hold on
%     loglog(eps_range(idx), vol1(idx), '.r')
%     hold off
%     pause
    P = polyfit(log(eps_range(idx)), log(vol1(idx)), 1);
    slopes(i) = P(1);
%     slopes(i,:) = (log(vol1(2:end))-log(vol1(1:end-1)))./(log(eps_range(2:end))-log(eps_range(1:end-1)));
    
end

figure;
plot(2*slopes)

% figure;
% plot(max(slopes,[],2))
% 
% figure;
% for i=1:nsamples
%     plot(slopes(i,:))
%     pause
% end


