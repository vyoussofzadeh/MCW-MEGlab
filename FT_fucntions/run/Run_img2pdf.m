d1 = rdir(source_img_dir);
for i=1:length(d1)
    imshow(imread(d1(i).name))
    set(gcf,'Color','w')
    export_fig(dest_pdf_dir,'-append')
end