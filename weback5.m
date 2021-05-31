% This code was not used in manuscript but may help 
% select a suitable choice of alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% f: Gray scale image 
%
% Output:
% maxAlpha: scalar value between -1.9 and 1.9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function maxAlpha = weback5(f)

%[M, N] = size(f);
n = 255;
tauAxis = 0:n;

alphaVec = -1.9:.05:1.9;

maxLength = 0;
maxAlpha = [];
for kk = 1:length(alphaVec)
    alphaIn = alphaVec(kk);
    % Compute scales
    Sf = localScale(f,alphaIn);
    % Compute averaged scales
    avgSf = mean(Sf,[2 3]);
    normSf = avgSf/max(avgSf(:));

    TV_Sf = sum(abs(diff(normSf)));
    
    if TV_Sf > maxLength 
        %update max length 
        maxAlpha = alphaIn;
        maxLength = TV_Sf;
    end
end