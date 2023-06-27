function d_in = do_readpow(cfg_main)

d_in  = cfg_main.d_in;

pow_LH = []; pow_RH = [];
for j=1:size(d_in.pow_sub,1)
    for i=1:size(d_in.pow_sub,2)
        pow_LH(j,i,:) = d_in.pow_sub(j,i).left;
        pow_RH(j,i,:) = d_in.pow_sub(j,i).right;
    end
end
d_in.pow_LH = pow_LH;
d_in.pow_RH = pow_RH;