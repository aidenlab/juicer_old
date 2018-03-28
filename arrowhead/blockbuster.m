function [result, scores, scores1] = blockbuster(observed, var_thresh, sign_thresh, list, list1)
% BLOCK_FINDER  Find blocks in input observed matrix.

    if (nargin < 4)
        list = []; 
        list1=[]; 
    end
    disp_figures=0;
    %calculate D upstream, directionality index upstream
    preferred_window=length(observed);
    gap=7;
    D_upstream=0.*observed;
    for i=1:length(observed)
        window=min(length(observed)-i-gap,i-gap);
        window=min(window,preferred_window);
        A=fliplr(observed(i,i-window:i-gap));
        B=observed(i,i+gap:i+window);
        
        preference=(A-B)./(A+B);
        D_upstream(i,i+gap:i+window)=preference;
    end
    
    % calculate triangles
    [Up, Up_Sign, Up_Sq, Lo, Lo_Sign, Lo_Sq] = dyp_triangles(D_upstream);
        
    B1=(-Up+Lo)/max(max(-Up+Lo))+(-Up_Sign+Lo_Sign)/max(max(-Up_Sign+Lo_Sign))-(Up_Sq-Up.^2+Lo_Sq-Lo.^2)/max(max(Up_Sq-Up.^2+Lo_Sq-Lo.^2));

    scores = zeros(size(list,1), 1);
    for i=1:size(list,1)
        scores(i) = max(max(B1(list(i,1):list(i,2), list(i,3):list(i,4))));
    end
    scores1 = zeros(size(list1,1), 1);
    for i=1:size(list1,1)
        scores1(i) = max(max(B1(list1(i,1):list1(i,2), list1(i,3):list1(i,4))));
    end
    
    B1(-Up_Sign < sign_thresh) = 0;
    B1(Lo_Sign  < sign_thresh)  = 0;
    Uvar = Up_Sq-Up.^2;
    Lvar = Lo_Sq-Lo.^2;
    if (var_thresh ~= 1000)
        B1(Uvar+Lvar > var_thresh)=0;
    end
    
    % find connected components of local max and average them
    
    CC1 = bwconncomp(B1>0);
    result = zeros(CC1.NumObjects, 7);
    
    for i=1:CC1.NumObjects
        [I, J] = ind2sub(CC1.ImageSize, CC1.PixelIdxList{i});
        [score, ind] = max(B1(CC1.PixelIdxList{i}));
        result(i, 1:2) = [I(ind),J(ind)];
        result(i,3) = score;
        result(i,4) = Uvar(I(ind),J(ind));
        result(i,5) = Lvar(I(ind),J(ind));
        result(i,6) = -Up_Sign(I(ind),J(ind));
        result(i,7) = Lo_Sign(I(ind),J(ind));
    end

    if (disp_figures)
        white2=cat(3, ones(size(D_upstream)), ones(size(D_upstream)),ones(size(D_upstream)));
        white=cat(3, zeros(size(D_upstream)), zeros(size(D_upstream)),ones(size(D_upstream)));
        red_color_map= makeColorMap([1 1 1],[1 0 0],80);
        mymap=zeros(size(observed));
        for k=1:size(result,1)
            row1=max(result(k,1)-1,1);
            row2=min(result(k,1)+1,length(observed));
            col1=max(result(k,2)-1,1);
            col2=min(result(k,2)+1,length(observed));
            
            mymap(row1:row2, row1:col2)=1;
            mymap(row1:col2, col1:col2)=1;
        end
        tmp=observed;  tmp(tmp>100)=100;
        h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(tmp)); hold off;
        set(h4, 'Position', [0 0 2000 2000]);
        colormap(red_color_map);
        set(h, 'AlphaData',mymap==0);
        h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(D_upstream)); hold off;
        set(h4, 'Position', [0 0 2000 2000]);
        h4=figure; h=imshow(white2, 'InitialMagnification', 300); axis('square'); %hold on; h=imshow(make_symmetric(B1)); hold off;
        set(h4, 'Position', [0 0 2000 2000]);
        set(h, 'AlphaData',B1~=0);
%         h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(Uvar)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
%         h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(Lvar)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
%          h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(-Up_Sign)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
%         h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(Lo_Sign)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
%          h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(-Up)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
%         h4=figure; imshow(white, 'InitialMagnification', 300); axis('square'); hold on; h=imagesc(make_symmetric(Lo)); hold off;
%         set(h4, 'Position', [0 0 2000 2000]);
    end
end
            
function Y = make_symmetric(X)
    Y = X;   
    for j=2:length(X)
        for i=1:j-1
            Y(j,i)=X(i,j);
        end
    end
end

function [Up, Up_Sign, Up_Sq, Lo, Lo_Sign, Lo_Sq] = dyp_triangles(M) 

  %calculate Bnew, the block score matrix. it's a combination of 3 matrices
    M(isnan(M))=0;
    window=length(M); % not using this because it messed things up

    % Matrices used as dynamic programming lookups.
    % "R" matrices are sums of the columns up to that point: R(1,5) is sum of
    % column 5 from diagonal (row 5) up to row 1
    % "U" matrices are sums of the rows up to the point: U(1,5) is sum of row 5
    % from diagonal (col 1) up to col 5
    % We want mean, mean of sign, and variance, so we are doing the sum then 
    % dividing by counts 
    Rsum=dyp_right(M, window);
    Rsign=dyp_right(sign(M), window);
    Rsq=dyp_right(M.*M, window);
    Rcount=dyp_right(ones(size(M)), length(M));

    Usum=dyp_upper(M,window);
    Usign=dyp_upper(sign(M), window);
    Usq=dyp_upper(M.*M,window);
    Ucount=dyp_upper(ones(size(M)), length(M));
    
    Up=0.*M;
    Up_Sign=0.*M;
    Up_Sq=0.*M;
    Up_Count=0.*M;
    for i=1:length(Up)
       for j=i+1:length(Up)
           bot = floor((j-i+1)/2);        
           % add half of column
           Up(i,j)= Up(i,j-1) + Rsum(i,j) - Rsum(i+bot, j);
           Up_Sign(i,j) = Up_Sign(i,j-1) + Rsign(i,j) - Rsign(i+bot, j);
           Up_Sq(i,j) = Up_Sq(i,j-1) + Rsq(i,j) - Rsq(i+bot, j);
           Up_Count(i,j)= Up_Count(i,j-1) + Rcount(i,j) - Rcount(i+bot, j);
       end
    end
   
    % Normalize
    Up_Count(Up_Count==0)=1;
    Up=Up./Up_Count;
    Up_Sign=Up_Sign./Up_Count;
    Up_Sq=Up_Sq./Up_Count;
    
    % Lower triangle
    Lo=0.*M;
    Lo_Sign=0.*M;
    Lo_Sq=0.*M;
    Lo_Count=0.*M;
    for a=1:length(Lo)
        for b=a+1:length(Lo)
            val=floor((b-a+1)/2);
            endpt = min(2*b-a,length(Lo));
            Lo_Count(a,b)=Lo_Count(a,b-1)+Ucount(b,endpt)-Rcount(a+val,b);
            Lo(a,b)=Lo(a,b-1)+Usum(b,endpt)-Rsum(a+val,b);
            Lo_Sign(a,b)=Lo_Sign(a,b-1)+Usign(b,endpt)-Rsign(a+val,b);
            Lo_Sq(a,b)=Lo_Sq(a,b-1)+Usq(b,endpt)-Rsq(a+val,b);
        end
    end
    Lo_Count(Lo_Count==0)=1;
    Lo=Lo./Lo_Count;
    Lo_Sign=Lo_Sign./Lo_Count;
    Lo_Sq=Lo_Sq./Lo_Count;
end
    
function M = flip_antidiagonal(N)
    M=N(end:-1:1,:);
    M=fliplr(M);
    M=M';
end

function U = dyp_upper(M, maxsize)
% DYP_UPPER Dynamic programming to calculate "upper" matrix
%   Initialize by setting the diagonal to the diagonal of original
%   Iterate down (for each row) and to the left.
    v=diag(M);
    U=diag(v);
    for i=1:length(U)
        endpoint=i+1+maxsize;
        if (endpoint > length(U))
            endpoint = length(U);
        end
        for j=i+1:endpoint
            U(i,j)=U(i,j-1)+M(i,j);
        end
    end
end

function R = dyp_right(M, maxsize)
% DYP_RIGHT Dynamic programming to calculate "right" matrix
%   Initialize by setting the diagonal to the diagonal of original
%   Iterate to the right and up.
    v=diag(M);
    R=diag(v);
    % j is column, i is row
    for j=2:length(R)
        endpoint=j-1-maxsize;
        if (endpoint < 1)
            endpoint = 1;
        end
        for i=j-1:-1:endpoint
            R(i,j)=M(i,j)+R(i+1,j);
        end
    end
end

function S = dyp_sum(M, super)
    v=diag(M);
    if (super > 1)
        S=zeros(size(M));
    else
        S=diag(v);
    end
    % d = distance from diagonal
    for d=super:length(S)-1
        % i = row, column is i+d
        % result is entry to left + entry below + orig value - diagonal
        % down (else we double count)
        for i=1:length(S)-d
            S(i,i+d)=S(i,i+d-1)+S(i+1,i+d)+M(i,i+d)-S(i+1,i+d-1);
        end
    end
end
    % flip over antidiagonal
    %observed2=flip_antidiagonal(observed);
    
    %calculate D downstream, directionality index downstream
%     D_downstream=0.*observed;
%     
%     
%     for i=1:length(observed2)
%         window=min(length(observed2)-i-gap,i-gap);
%         window=min(window,preferred_window);
%         A=fliplr(observed2(i,i-window:i-gap));
%         B=observed2(i,i+gap:i+window);
%         
%         preference=(A-B)./(A+B);
%         D_downstream(i,i+gap:i+window)=preference;
%     end   
% 
%     [Right, Right_Sign, Right_Sq, Left, Left_Sign, Left_Sq] = dyp_triangles(D_downstream);
%     
%     % flip back over antidiagonal
%     Right = flip_antidiagonal(Right);
%     Right_Sign = flip_antidiagonal(Right_Sign);
%     Right_Sq = flip_antidiagonal(Right_Sq);
%     Left = flip_antidiagonal(Left);
%     Left_Sign = flip_antidiagonal(Left_Sign);
%     Left_Sq = flip_antidiagonal(Left_Sq);
  
%     gt=[370 424
% 368 396
% 425 450
% 450 462
% 425 462
% 465 468
% 468 496
% 491 496
% 496 503
% 503 516
% 540 672
% 673 707
% 707 713
% 713 720
% 708 731
% 738 763
% 768 1107
% 768 940
% 940 981
% 981 1104
% 1112 1207
% ];
    
    
%     B2=(-Right+Left)+(-Right_Sign+Left_Sign)-(Right_Sq-Right.^2+Left_Sq-Left.^2);
%     B2(Left_Sign < sign_thresh) = 0;
%     B2(Right_Sign > -sign_thresh) = 0;
%     Rvar = Right_Sq-Right.^2;
%     Levar = Left_Sq-Left.^2;  
%     B2(Rvar > var_thresh)=0;
%     B2(Levar > var_thresh)=0;
%     for i=1:CC2.NumObjects
%         [I, J] = ind2sub(CC2.ImageSize, CC2.PixelIdxList{i});
%         [score, ind] = max(B2(CC2.PixelIdxList{i}));
%         result(i+CC1.NumObjects, 1:2) = [I(ind),J(ind)];
%         result(i+CC1.NumObjects,3) = score;
%         result(i+CC1.NumObjects,4) = Utotal(I(ind),J(ind));
%         result(i+CC1.NumObjects,5) = Ltotal(I(ind),J(ind));
%         result(i+CC1.NumObjects,6) =  Uvar(I(ind),J(ind));
%         result(i+CC1.NumObjects,7) = Lvar(I(ind),J(ind));
%         result(i+CC1.NumObjects,8) = Letotal(I(ind),J(ind));
%         result(i+CC1.NumObjects,9) = Rtotal(I(ind),J(ind));
%         result(i+CC1.NumObjects,10) = Rvar(I(ind),J(ind));
%         result(i+CC1.NumObjects,11) = Levar(I(ind),J(ind)) ; 
%     end
%     count =0;
%     for i=1:length(gt) 
%         a=gt(i,1); b=gt(i,2);
%         %disp([ B1(a,b) B2(a,b) Utotal(a,b) Ltotal(a,b) Uvar(a,b) Lvar(a,b) Letotal(a,b) Rtotal(a,b) Rvar(a,b) Levar(a,b)]);
%         if (B1(a,b)>0)
%             count=count+1;
%         end
%     end
%     
%     precision = count/CC1.NumObjects;
%     recall = count/length(gt);
%     disp([2*(precision*recall/(precision+recall)),precision, recall]);
%     


