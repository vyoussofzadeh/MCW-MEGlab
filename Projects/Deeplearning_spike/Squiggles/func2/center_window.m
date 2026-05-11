% helper
function y = center_window(x, W)
  g = std(x,[],1); [~,i] = max(g); T = size(x,2);
  L = max(1, i - floor(W/2)); R = min(T, L + W - 1);
  y = x(:, L:R);
end
