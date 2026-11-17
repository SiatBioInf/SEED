d = read.csv('D:/科研/CTC知识图谱/v2/data/processed/gene_cancer_full_stats.csv')
epcam = d[d$gene == 'EPCAM', c('cancer','mean_percentile','pct_delta','detection_rate','confidence_tier')]
vim = d[d$gene == 'VIM', c('cancer','mean_percentile','pct_delta','detection_rate','confidence_tier')]
names(epcam)[2:5] = c('EPCAM_pct','EPCAM_delta','EPCAM_det','EPCAM_tier')
names(vim)[2:5] = c('VIM_pct','VIM_delta','VIM_det','VIM_tier')
m = merge(epcam, vim, by='cancer')
m = m[order(m$EPCAM_pct, decreasing=TRUE), ]
cat(sprintf('%-6s | %9s %10s %9s | %7s %9s %8s | %s\n',
    'Cancer','EPCAM%','EPCAM_d','EPCAM_det','VIM%','VIM_d','VIM_det','Phenotype'))
cat(strrep('-', 90), '\n')
for(i in 1:nrow(m)) {
  r = m[i,]
  if(r$EPCAM_pct < 0.6 && r$VIM_pct > 0.7) {
    phenotype = 'MESENCHYMAL (EMT)'
  } else if(r$EPCAM_pct > 0.7 && r$VIM_pct < 0.7) {
    phenotype = 'EPITHELIAL'
  } else if(r$EPCAM_pct > 0.7 && r$VIM_pct > 0.7) {
    phenotype = 'HYBRID'
  } else {
    phenotype = 'intermediate'
  }
  cat(sprintf('%-6s | %9.3f %10.3f %9.2f | %7.3f %9.3f %8.2f | %s\n',
      r$cancer, r$EPCAM_pct, r$EPCAM_delta, r$EPCAM_det,
      r$VIM_pct, r$VIM_delta, r$VIM_det, phenotype))
}
