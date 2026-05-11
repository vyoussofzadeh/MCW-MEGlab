function Xd = decimate_time_2x_4d(X) % X: [C x T x Cin x N]
Xd = (X(:,1:2:end,:,:)+X(:,2:2:end,:,:))/2; % avg-pool stride 2
end
