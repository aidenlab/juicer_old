function [final, scores] = run_blockbuster(M,str,list, control)

  if (nargin < 3)
      list=[];
      control=[];
  end
  % get large number of blocks
  [result1, scores, control_scores]=call_blockbuster(M, 1000, 0.4, list, control);
  p_sgn=0.4;
  while (size(result1,1)==0 && p_sgn > 0)
    p_sgn = p_sgn-0.1;
    [result1, scores, control_scores]=call_blockbuster(M, 1000, p_sgn, list, control);
  end
 
  % high variance threshold, fewer blocks, high confidence
  result2 = call_blockbuster(M, 0.2, 0.5);
  
  resdiff = setdiff(result1,result2, 'rows');
  % only interested in big blocks.
  result1=resdiff(resdiff(:,2)-resdiff(:,1)>60, :);
  % result2 is in.  then take only bigger blocks in "block deserts"
  for i=1:size(result1,1)
      conflicts=false;
      for j=1:size(result2,1)
          % if the beginning of an existing block is within this interval,
          % or the end is within this interval, this is a conflicting block
          if ((result2(j,1) >= result1(i,1) && result2(j,1) <= result1(i,2)) || (result2(j,2) >= result1(i,1) && result2(j,2) <= result1(i,2)))
              conflicts=true;
          end
      end
      if (~conflicts)
          result2(end+1,:) = result1(i,:);
      end
  end
  
  result=result2;
  
  fprintf('\n');
  scores=merge_scores(scores);
  control_scores=merge_scores(control_scores);
  
  if (size(result,1)==0)
      final = [];
      return;
  end
     
  s = struct;
  s.minx=result(1,1); s.maxx=result(1,1);
  s.miny=result(1,2); s.maxy=result(1,2);
  s.value = result(1,:);
  C{1} = s;
  for i=2:size(result,1)
    added=false;

    for j=1:length(C)
        if ((abs(C{j}.minx - result(i,1)) < 5 || abs(C{j}.maxx - result(i,1)) < 5) && (abs(C{j}.miny - result(i,2)) < 5 || abs(C{j}.maxy - result(i,2)) < 5))
            % add to this bin
            added=true;
            C{j}.value(end+1,:) =result(i,:);
            if (result(i,1) < C{j}.minx); C{j}.minx = result(i,1); end
            if (result(i,1) > C{j}.maxx); C{j}.maxx = result(i,1); end
            if (result(i,2) < C{j}.miny); C{j}.miny = result(i,2); end
            if (result(i,2) > C{j}.maxy); C{j}.maxy = result(i,2); end
            break; % don't look anymore in C for bins
        end
    end
    if (added ~= true)
        s=struct;
        s.value=result(i,:);
        s.minx=result(i,1);
        s.maxx=result(i,1);
        s.miny=result(i,2);
        s.maxy=result(i,2);
        C{end+1}=s;
    end
  end
    final = zeros(length(C),4);
    for i=1:length(C)
        final(i,1)=C{i}.maxx;
        final(i,2)=C{i}.maxy;
        final(i,3)=mean(C{i}.value(:,3));
        final(i,4)=mean(C{i}.value(:,4));
        final(i,5)=mean(C{i}.value(:,5));
        final(i,6)=mean(C{i}.value(:,6));
        final(i,7)=mean(C{i}.value(:,7));
    end

  result=final;
  clear C;
  s = struct;
  s.minx=result(1,1); s.maxx=result(1,1);
  s.miny=result(1,2); s.maxy=result(1,2);
  s.value = result(1,:);
  C{1} = s;
  for i=2:size(result,1)
    added=false;

    for j=1:length(C)
        if ((abs(C{j}.minx - result(i,1)) < 10 || abs(C{j}.maxx - result(i,1)) < 10) && (abs(C{j}.miny - result(i,2)) < 10 || abs(C{j}.maxy - result(i,2)) < 10))
            % add to this bin
            added=true;
            C{j}.value(end+1,:) =result(i,:);
            if (result(i,1) < C{j}.minx); C{j}.minx = result(i,1); end
            if (result(i,1) > C{j}.maxx); C{j}.maxx = result(i,1); end
            if (result(i,2) < C{j}.miny); C{j}.miny = result(i,2); end
            if (result(i,2) > C{j}.maxy); C{j}.maxy = result(i,2); end
            break; % don't look anymore in C for bins
        end
    end
    if (added ~= true)
        s=struct;
        s.value=result(i,:);
        s.minx=result(i,1);
        s.maxx=result(i,1);
        s.miny=result(i,2);
        s.maxy=result(i,2);
        C{end+1}=s;
    end
  end
final = zeros(length(C),4);
for i=1:length(C)
    final(i,1)=C{i}.maxx;
    final(i,2)=C{i}.maxy;
    final(i,3)=mean(C{i}.value(:,3));
    final(i,4)=mean(C{i}.value(:,4));
    final(i,5)=mean(C{i}.value(:,5));
    final(i,6)=mean(C{i}.value(:,6));
    final(i,7)=mean(C{i}.value(:,7));
end

[~,ind]=sort(final(:,4)+final(:,5), 'descend');
 
    
final = final(ind,:);
dlmwrite(str,final);
if (~isempty(list))
    dlmwrite([str '_scores'], scores);
    dlmwrite([str '_control_scores'], control_scores);
end

function [result, scores, control_scores] = call_blockbuster(M, p_var, p_sign, list, control)
  if (nargin < 4)
      list=[];
  end
  result=[];
  scores=[];
  control_scores=[];

  endpt = size(M,1);
  for limstart=1:1000:endpt
    limend = limstart+2000;
    if (limend > endpt)
        limend = endpt;
    end
    if (~isempty(list))
        lst = list(list(:,1) >= limstart & list(:,1) <= limend & list(:,2) >= limstart & list(:,2) <= limend & list(:,3) >= limstart & list(:,3) <= limend & list(:,4) >= limstart & list(:,4) <= limend, :);
        lst1 = control(control(:,1) >= limstart & control(:,1) <= limend & control(:,2) >= limstart & control(:,2) <= limend & control(:,3) >= limstart & control(:,3) <= limend & control(:,4) >= limstart & control(:,4) <= limend, :);
    else
        lst=[]; 
        lst1=[];
    end

    [M1, sr, sr1]=blockbuster(full(M(limstart:limend,limstart:limend)),p_var,p_sign, lst-limstart+1, lst1-limstart+1);

    if (~isempty(M1))
        M1(:, 1:2)=M1(:,1:2)+limstart-1;
    end
    result=[result; M1];
    scores=[scores; lst, sr];
    control_scores=[control_scores; lst1, sr1];
    
    fprintf('.');
end

function res = merge_scores(scores)
  scores=sortrows(scores);
  for i=1:size(scores,1)
      if (isnan(scores(i,5)))
          scores(i,5)=-10000;
      end
      for j=i+1:size(scores,1)
          if (scores(i,1:4) == scores(j,1:4))
              scores(i,5)=nanmax(scores(i,5), scores(j,5));
              scores(j,5)=scores(i,5);
          end
      end
  end
  res=unique(scores, 'rows');
res(res==-10000)=NaN;

