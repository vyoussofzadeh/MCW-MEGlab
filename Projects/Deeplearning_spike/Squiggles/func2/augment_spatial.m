function [I_aug, Y_aug] = augment_spatial(Icells, Y, rot90on, doFlip)
I_aug = {}; Y_aug = categorical();
for i=1:numel(Icells)
    im = Icells{i}; cls=Y(i);
    I_aug{end+1,1}=im; Y_aug(end+1,1)=cls;
    if rot90on
        I_aug{end+1,1}=rot90(im,1); Y_aug(end+1,1)=cls;
        I_aug{end+1,1}=rot90(im,2); Y_aug(end+1,1)=cls;
        I_aug{end+1,1}=rot90(im,3); Y_aug(end+1,1)=cls;
    end
    if doFlip
        I_aug{end+1,1}=-im; Y_aug(end+1,1)=cls;
    end
end
end