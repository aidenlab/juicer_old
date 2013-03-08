function plotstats(A,B,D,x,namestr)
  h = figure;
  g = semilogy(A./sum(A));
  set(g, 'LineWidth', 1.5);
  set(gca, 'FontSize', 16);
  xlabel('Distance (bp)');
  ylabel('Fraction of Reads (log)');
  title('Distance from closest RE site');
  set(gcf, 'PaperPositionMode', 'auto');
  str = sprintf('%s_RE.jpg', namestr);
  print(h, '-djpeg100', str);
  h = figure;
  for i=1:size(B,1)
    mapqthresh(i,1) = sum(B(i:end,1));
    mapqthresh(i,2) = sum(B(i:end,2));
    mapqthresh(i,3) = sum(B(i:end,3));
  end
    g = plot(0:200, mapqthresh(1:201,1), '-b', 0:200, mapqthresh(1:201,2), '-r', 0:200, mapqthresh(1:201,3), '-g');
%  g = bar(0:200,B);
  set(g, 'LineWidth', 1.5);
  xlim([0 200]);
  set(gca, 'FontSize', 16);
  xlabel('MapQ threshold');
  ylabel('Count');
  title('MapQ threshold count');
  legend('All', 'Intra', 'Inter');
  str = sprintf('%s_MAPQ.jpg', namestr);
  set(gcf, 'PaperPositionMode', 'auto');
  print(h, '-djpeg100', str);
  h = figure;
  g = loglog(x,D(:,1),'-r',x,D(:,2),'-g',x,D(:,3),'-b',x,D(:,4),'-c');
  set(g, 'LineWidth', 1.5);
  set(gca, 'FontSize', 16);
  xlim([10, 10^10]);
  legend('Inner', 'Outer', 'Right', 'Left', 'Location', 'Best');
  ylabel('Number of Binned Reads (log)');
  xlabel('Distance (log)');
  title('Types of reads vs distance');
  str = sprintf('%s_readType.jpg', namestr);
  set(gcf, 'PaperPositionMode', 'auto');
  print(h, '-djpeg100', str);
  h = figure;
  E = sum(D,2);
  F = cumsum(E);
  g = semilogx(x,F);
  set(g, 'LineWidth', 1.5);
  set(gca, 'FontSize', 12);
  hold on;
  g = semilogx([x(37),x(37)],[min(F),max(F)],'-r');
  set(g, 'LineWidth', 1.5);
  g = semilogx([x(56),x(56)],[min(F),max(F)],'-g');
  set(g, 'LineWidth', 1.5);
	str1(1)={['\fontsize{12}' sprintf('%d less than', F(37))]};
  str1(2)={sprintf('%d apart', x(37))};
text(x(37),max(F),str1, 'HorizontalAlignment', 'right')
	str2(1)={['\fontsize{12}' sprintf('%d less than', F(56))]};
  str2(2)={sprintf('%d apart', x(56))};
text(x(56),max(F),str2)
  ylabel('Cumulative Sum of Binned Reads (log)');
  xlabel('Distance (log)');
  title('Intra reads vs distance');
  str = sprintf('%s_intraCount.jpg', namestr);
  set(gcf, 'PaperPositionMode', 'auto');
  print(h, '-djpeg100', str);
end
