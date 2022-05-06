% save
if flag.savepdf
    set(gcf,'Color','w')
    export_fig(['subject_pdfs/' [subj,'_', run] '.pdf'],'-append')
end


