load SI_allelicCN_1Mb.mat ;
COV=AllelicCN_1Mb;
SEG=AllelicSeg;
Samples=SampleInfo.SampleNames;
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Seg=SEG{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
		chrlist=unique(Cov.chr);
		for chri=1:length(chrlist)
			chrlabel=chrlist{chri};
			idx=strcmp(Cov.chr,chrlabel);
			plotPos=Cov.pos(idx);
			plotAlleleA=Cov.alleleA(idx);
			plotAlleleB=Cov.alleleB(idx);
			ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
			hold on,
			seg=Seg.B(strcmp(Seg.B.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;
			for i=1:length(seg)-1
				plot([seg.Start(i),seg.End(i)],seg.AvgDepth(i)*[1,1],'r','LineWidth',1.5);
				if (seg.End(i)==seg.Start(i+1))
					plot(seg.End(i)*[1,1],[seg.AvgDepth(i),seg.AvgDepth(i+1)],'r','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'r','LineWidth',1.5);
			
			seg=Seg.A(strcmp(Seg.A.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;

			for i=1:length(seg)-1
				plot([seg.Start(i),seg.End(i)],seg.AvgDepth(i)*[1,1],'b','LineWidth',1.5);
				if (seg.End(i)==seg.Start(i+1))
					plot(seg.End(i)*[1,1],[seg.AvgDepth(i),seg.AvgDepth(i+1)],'b','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'b','LineWidth',1.5);

			%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
			print(gcf,'-dpdf',[Samples{i} '_' chrlabel '_1Mb.pdf']);
			close all
		end
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end
		
