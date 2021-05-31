function f = yellowCircles(~)

N=256;  % some size of grid
if mod(N,2) % odd vs even matrix sizes
        [x,y] = meshgrid(-(N-1)/2:(N-1)/2);
else
        [x,y] = meshgrid(-N/2+1:N/2);
end
% assuming the circle are always in the center, but we can modify this if needed 
x0=0; y0=0; 

% say we want a ring between `r1` and `r2`
Circ = @(r1,r2) (x-x0).^2+(y-y0).^2<=r2^2 & ...
             (x-x0).^2+(y-y0).^2>=r1^2;
         
% imagesc(Circ(30,40)+Circ(50,60))

f = Circ(10,20)+Circ(50,60)+...
    Circ(90,100); 