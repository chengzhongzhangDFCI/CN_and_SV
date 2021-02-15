%% plot raw sequence coverage, segmented total coverage, and segmental average allelic coverage
load hg38.mat
for i=1:length(SampleNames)
	chrlabel='chr4';
	maxLength=6e7;
	maxLength=max(bins.End);
	maxX=log10(double(maxLength)/1e6);
	XScale=maxX*500;
    pos=double(bins.End)/1e6;
	cov=chr4Cov{i}.gc_cov;
	maxY=ceil(max(cov));
	YScale=maxY*50;
	figure,set(gcf, 'Unit', 'pixels', 'Position', [15,0,XScale+200 YScale+75]),hold on;
	plot(pos,cov,'o','MarkerFaceColor',[0.75,0.75,0.75],'MarkerEdgeColor',[0.75,0.75,0.75],'MarkerSize',3);	
	seg=chr4Seg{i};
	seg.Start=double(seg.Start)/1e6;seg.End=double(seg.End)/1e6;
	plot([seg.Start(1),seg.End(1)],[1,1]*seg.seg_total(1),'k-','LineWidth',3);
	plot([seg.Start(1),seg.End(1)],[1,1]*seg.seg_alleleA(1),'r-','LineWidth',3);
	plot([seg.Start(1),seg.End(1)],[1,1]*seg.seg_alleleB(1),'b-','LineWidth',3);
	for j=2:length(seg)
		plot([seg.Start(j),seg.Start(j)],[0 maxY-1],'k--');
		plot([seg.Start(j),seg.End(j)],[1,1]*seg.seg_total(j),'k-','LineWidth',3);
		plot([seg.Start(j),seg.End(j)],[1,1]*seg.seg_alleleA(j),'r-','LineWidth',3);
		plot([seg.Start(j),seg.End(j)],[1,1]*seg.seg_alleleB(j),'b-','LineWidth',3);
	end

	gBands=GiemsaBands(strcmp(GiemsaBands.Chr,'chr4'),:);
	gBands.Start=double(gBands.Start)/1e6;gBands.End=double(gBands.End)/1e6;
	gBandColor=gBands.Stain;
	gBandColor(strcmp(gBandColor,'acen'))={[0.9,0,0]};
	gBandColor(strcmp(gBandColor,'gpos25'))={[0.75,0.75,0.75]};
	gBandColor(strcmp(gBandColor,'gpos50'))={[0.5,0.5,0.5]};
	gBandColor(strcmp(gBandColor,'gpos75'))={[0.25,0.25,0.25]};
	gBandColor(strcmp(gBandColor,'gpos100'))={[0,0,0]};
	gBandColor(strcmp(gBandColor,'stalk'))={[0,0.5,0]};
	gBandColor(strcmp(gBandColor,'gvar'))={[0.5,0.25,0]};
	gBandColor(strcmp(gBandColor,'gneg'))={[1,1,1]};
	rectangle('Position',[gBands.Start(1),-0.95,gBands.End(end)-gBands.Start(1),0.6],'LineWidth',1.5);
	for gi=1:length(gBands)
		rectangle('Position',[gBands.Start(gi),-0.95,gBands.End(gi)-gBands.Start(gi),0.6],'FaceColor',gBandColor{gi},'EdgeColor','none');
	end
	box off,set(gca,'TickDir','out','FontSize',16,'Ylim',[-0.95 maxY],'Xlim',[-0.02*pos(end),ceil(maxLength/1e7)*10],'XMinorTick','On');
	set(gca,'Unit','pixels','Position',[100 50 XScale YScale]);
	title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
	set(gcf, 'Unit', 'inches');
	Position=get(gcf,'Position');
	set(gcf, 'PaperUnits', 'inches', 'PaperSize', Position(3:4), 'PaperPosition', [0,0,Position(3:4)]);
	print(gcf,'-dpdf',[SampleNames{i} '_' chrlabel '_10kb.pdf']);
	close all
end
