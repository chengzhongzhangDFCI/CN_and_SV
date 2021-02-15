load SI_allelicCN_1Mb.mat;
COV=AllelicCN_1Mb;
load SampleSV.mat
RA=SampleSV_Depth8;
Samples=SampleInfo.SampleNames;
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
		SV=RA{i};
		SV=SV(SV.TotalCount-SV.SplitCount>0 & (SV.SplitCount>0 | strcmp(SV.chr1,SV.chr2)),:);
		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
		export2circos([lrbkps;itbkps],Cov,[Samples{i}]);	
		export2circos([lrbkps;itbkps],[],[Samples{i}]);		
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end

