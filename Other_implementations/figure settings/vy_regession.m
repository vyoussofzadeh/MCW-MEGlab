tbl = table(X(:,idx(1)),y1');

mdl = fitlm(tbl,'linear');
plot(mdl)