load SI_allelicCN_1Mb.mat;
COV=AllelicCN_1Mb;
load SampleSVs.mat;
load SampleInfo.mat;
Samples=SampleInfo.SampleNames;
cd CIRCOS
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
		SV=SampleSVs{i};
		SV=SV(SV.suppSamples==1 & SV.NCount==0,:);
		SV=SV(SV.TotalCount>2,:);
		SV=SV((SV.TotalCount-SV.SplitCount>0 & SV.SplitCount>0) | strcmp(SV.chr1,SV.chr2),:);	
		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
		export2circos([lrbkps;itbkps],Cov,[Samples{i}]);	
		%export2circos([lrbkps;itbkps],[],[Samples{i}]);		
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end
