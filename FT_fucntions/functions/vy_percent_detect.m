function l = vy_percent_detect(xin)

xin(isnan(xin)) = [];
mxin = mean(xin);
l = length(find(xin >= mxin))/length(xin);

end