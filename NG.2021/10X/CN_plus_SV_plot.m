load SampleInfo.mat
load SI_allelicCN_1Mb.mat ;
COV=AllelicCN_1Mb;
load AllelicCNSeg.mat
SEG=AllelicSeg;
SEG2=AllelicSeg2;
load SampleSVs.mat

RA=SampleSVs;
Samples=SampleInfo.SampleNames;
if ~exist('HiThreshold_1Mb')
	mkdir('HiThreshold_1Mb');
end
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
		Seg=SEG{i};
		Seg2=SEG2{i};
		SV=RA{i};
		SV=SV(SV.TotalCount>2,:);
		SV=SV((SV.TotalCount-SV.SplitCount>0 & SV.SplitCount>0) | strcmp(SV.chr1,SV.chr2),:);
		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
		lrbkps.pos1=double(lrbkps.pos1)/1e6;lrbkps.pos2=double(lrbkps.pos2)/1e6;
		itbkps.pos1=double(itbkps.pos1)/1e6;itbkps.pos2=double(itbkps.pos2)/1e6;
		chrlist=unique(Cov.chr);
		for chri=1:length(chrlist)
			chrlabel=chrlist{chri};
			idx=strcmp(Cov.chr,chrlabel);
			plotPos=Cov.pos(idx);
			plotAlleleA=Cov.alleleA(idx);
			plotAlleleB=Cov.alleleB(idx);
			ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
			hold on,
			% diploid cell segmentation
			seg=Seg.B(strcmp(Seg.B.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;
			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1],'r','LineWidth',1.5);
				if (seg.End(si)==seg.Start(si+1))
					plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)],'r','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'r','LineWidth',1.5);
			
			seg=Seg.A(strcmp(Seg.A.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;

			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1],'b','LineWidth',1.5);
				if (seg.End(si)==seg.Start(si+1))
					plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)],'b','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'b','LineWidth',1.5);
			
			% G2 cell segmentation (tetraploid)
			seg=Seg2.B(strcmp(Seg2.B.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;
			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1]/2,'r--','LineWidth',1.5);
				%if (seg.End(si)==seg.Start(si+1))
				%	plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)]/2,'r','LineWidth',0.5);
				%end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1]/2,'r--','LineWidth',1.5);
			
			seg=Seg2.A(strcmp(Seg2.A.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;

			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1]/2,'b--','LineWidth',1.5);
				%if (seg.End(si)==seg.Start(si+1))
			%		plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)]/2,'b','LineWidth',0.5);
			%	end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1]/2,'b--','LineWidth',1.5);

			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1~=lrbkps.str2,:),4,1,'k');
			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1==lrbkps.str2,:),4,-1,'k');
			bkps1=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==-1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==-1)];
			bkps2=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==1)];
			for ki=1:length(bkps1)
				plot(bkps1(ki)*[1,1],[4,5],'m');
			end
			for ki=1:length(bkps2)
				plot(bkps2(ki)*[1,1],[4,3],'m');
			end
			%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
			print(gcf,'-dpdf',['HiThreshold_1Mb/' Samples{i} '_' chrlabel '_1Mb.pdf']);
			close all
		end
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end

if ~exist('LowThreshold_1Mb')
	mkdir('LowThreshold_1Mb');
end
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
		Seg=SEG{i};
		Seg2=SEG2{i};
		SV=RA{i};
		SV=SV((SV.TotalCount-SV.SplitCount>0 & SV.SplitCount>0 & SV.TotalCount>2) | strcmp(SV.chr1,SV.chr2),:);
		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
		lrbkps.pos1=double(lrbkps.pos1)/1e6;lrbkps.pos2=double(lrbkps.pos2)/1e6;
		itbkps.pos1=double(itbkps.pos1)/1e6;itbkps.pos2=double(itbkps.pos2)/1e6;
		chrlist=unique(Cov.chr);
		for chri=1:length(chrlist)
			chrlabel=chrlist{chri};
			idx=strcmp(Cov.chr,chrlabel);
			plotPos=Cov.pos(idx);
			plotAlleleA=Cov.alleleA(idx);
			plotAlleleB=Cov.alleleB(idx);
			ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
			hold on,
			% diploid cell segmentation
			seg=Seg.B(strcmp(Seg.B.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;
			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1],'r','LineWidth',1.5);
				if (seg.End(si)==seg.Start(si+1))
					plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)],'r','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'r','LineWidth',1.5);
			
			seg=Seg.A(strcmp(Seg.A.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;

			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1],'b','LineWidth',1.5);
				if (seg.End(si)==seg.Start(si+1))
					plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)],'b','LineWidth',0.5);
				end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1],'b','LineWidth',1.5);
			
			% G2 cell segmentation (tetraploid)
			seg=Seg2.B(strcmp(Seg2.B.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;
			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1]/2,'r--','LineWidth',1.5);
				%if (seg.End(si)==seg.Start(si+1))
				%	plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)]/2,'r','LineWidth',0.5);
				%end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1]/2,'r--','LineWidth',1.5);
			
			seg=Seg2.A(strcmp(Seg2.A.Chr,chrlabel),:);
			seg.Start=double(seg.Start)/1e6;
			seg.End=double(seg.End)/1e6;

			for si=1:length(seg)-1
				plot([seg.Start(si),seg.End(si)],seg.AvgDepth(si)*[1,1]/2,'b--','LineWidth',1.5);
				%if (seg.End(si)==seg.Start(si+1))
			%		plot(seg.End(si)*[1,1],[seg.AvgDepth(si),seg.AvgDepth(si+1)]/2,'b','LineWidth',0.5);
			%	end
			end
			plot([seg.Start(end),seg.End(end)],seg.AvgDepth(end)*[1,1]/2,'b--','LineWidth',1.5);

			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1~=lrbkps.str2,:),4,1,'k');
			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1==lrbkps.str2,:),4,-1,'k');
			bkps1=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==-1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==-1)];
			bkps2=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==1)];
			for ki=1:length(bkps1)
				plot(bkps1(ki)*[1,1],[4,5],'m');
			end
			for ki=1:length(bkps2)
				plot(bkps2(ki)*[1,1],[4,3],'m');
			end
			%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
			print(gcf,'-dpdf',['LowThreshold_1Mb/' Samples{i} '_' chrlabel '_1Mb.pdf']);
			close all
		end
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end


