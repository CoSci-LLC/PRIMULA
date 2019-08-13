% ProbSLOPE (Cislaghi 2016)

function [ProbSLOPE] = PrSlope(SLOPE,w,l)

meanA = NaN(size(SLOPE));
varA  = NaN(size(SLOPE));
nr = size(SLOPE,1); 
nc = size(SLOPE,2);

addw = floor((w/5-1)/2); 
addl = floor((l/5-1)/2);

for n = 1:numel(SLOPE)
    if ~isnan(SLOPE(n))
        [~,left] = ind2sub(size(SLOPE),n-nr*addw);
        if left <= 0 
            left = 1;
        end
        [~,right] = ind2sub(size(SLOPE),n+nr*addw);
        if right > size(SLOPE,2)
            right = size(SLOPE,2);
        end
        [row,~] = ind2sub(size(SLOPE),n);
        up = row - addl;
        down = row + addl;
        
        if up <= 0
            up = 1;
        end        
        
        if down > size(SLOPE,1)
            down = size(SLOPE,1);
        end
        B = SLOPE(up:down, left:right);
        meanA(n) = nanmean(B(:));
        varA(n)  = nanvar(B(:));
    end  
end

ProbSLOPE = arrayfun(@(x,y) normrnd(x,y),meanA,varA);
