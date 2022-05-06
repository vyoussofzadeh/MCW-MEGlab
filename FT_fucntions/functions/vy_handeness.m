function [HI] = vy_handeness(datain, thre)


HI = [];
HI.l = 0; HI.b = 0; HI.r = 0;
for i=1:length(datain)
    if datain(i) >= thre
        HI.l =  HI.l+1;
        hidx(i) = 1;
    elseif datain(i) < thre && datain(i) > -thre
        HI.b =  HI.b+1;
        hidx(i) = 2;
    elseif datain(i) <= -thre
        HI.r =  HI.r+1;
        hidx(i) = 3;
    end
end

HI.result = [datain;hidx];

HI.lp = HI.l/i*100;
HI.bp = HI.b/i*100;
HI.rp = HI.r/i*100;


