function GenomeWideCovPlot(Cov)
ColorOptions.option='ro';
ScatterPlot(Cov(:,1:2),Cov.alleleB,'Yscale',4,'Scalebar','Off','Labels','Off','ColorOptions',ColorOptions,'Normalize','Off');
ColorOptions.option='bo';
ScatterPlot(Cov(:,1:2),Cov.alleleA,'Yscale',4,'ColorOptions',ColorOptions,'Normalize','Off');
%plot_alleleCN_logTransformed_v2(loci,refCov(:,i),altCov(:,i),Haplotype,binSize);
 
