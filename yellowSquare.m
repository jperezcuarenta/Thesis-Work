function f = yellowSquare(~)
%%%%%%%%%%%%%
% Build Image
%%%%%%%%%%%%%


% top left
f1 = false(128,128);

% top right 
% % % Gives yellow scales
f21 = 1:128;
f22 = f21+1;
f21 = mod(f21,2);
f22 = mod(f22,2);
f2 = repmat([f21;f22],64,1);

% % does not give yellow scales
% z = binarySquare(32);
% f2 = repmat(z,4,4);

% bottom right
z = binarySquare(64);
f3 = repmat(z,2,2);

% bottom left
f4 = true(128,128);

% build 256 image
f = [f1 f2; f4 f3];
end
