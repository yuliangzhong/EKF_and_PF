% N = 5;
% kappa = 2 * 0.7 * rand(1,N) - 0.7;
% posx = (3-kappa) .* rand(1,N) + kappa;
% posy = 3 * rand(1,N);

pos = [1,nan,3,nan,5,11,12,13,nan,14;
       6,nan,8,nan,10,15,16,17,nan,18];
k = ones(1,10);

posx = pos(1,:);

ind = ~isnan(pos(1,:));

pp = posx(ind);
% 
% len = length(P);
% PP = reshape(P,[2,len/2]);
% k = k(~isnan(pos));











% beta = [0, 0.05,0.1,0.15,0.2,0.2,0.15,0.1,0.05];
% N = 9;
% Sigma = cumsum(beta);
% ind = zeros(1,N);
% for i = 1:N
%     r = rand();
%     j = N;
%     while r<= Sigma(j-1)
%         j = j-1;
%     end
%     ind(i) = j;
% end