clc;
close all;
clear all;

disp('Fitting a linear model to the termperature-curvature data in');
disp('Ranatunga 1982 (Table 1)')
fprintf('\n')
disp('Ranatunga KW. Temperature‐dependence of shortening velocity and rate of ');
disp('isometric tension development in rat skeletal muscle. The Journal of ');
disp('Physiology. 1982 Aug 1;329(1):465-83');
fprintf('\n')

t=[35, 30, 25, 20]';
c=[0.212, 0.18, 0.157, 0.137]';

% c0 + c1*t = c
% A*x = b 
% [35, 1](c1) = 0.212 
% [30, 1](c0) = 0.18 
% [25, 1]     = 0.157 
% [20, 1]     = 0.137

A = [t , ones(size(t))];
b = c;

% Ax = b
% A'Ax = A'b
% x = (A'A)\(A'b)

x = (A'*A)\(A'*b);

c12 = [12,1]*x;

disp('Expected curvature value at 12 C in whole soleus muscle');
disp(c12);
fprintf('\n')

disp('Solving for fv at vceMax*0.5');

vmax = 1.02;
% vmax for YF
%
% Degens H, Yu F, Li X, Larsson L. Effects of age and gender on shortening 
% velocity and myosin isoforms in single rat muscle fibres. Acta physiologica 
% scandinavica. 1998 May;163(1):33-40.

Po = 1;
c = c12;
b  = c*vmax;

v = [0.01:0.01:1].*vmax;
fv = ((1+c)*b - c.*(v+b)) ./ (v+b);

vHalf = 0.5*vmax;
fvHalf = ((1+c)*b - c.*(vHalf+b)) ./ (vHalf+b);

disp('fv at 0.5*vceMax');
disp(fvHalf);

plot(v,fv);
hold on
plot(vHalf,fvHalf,'o','MarkerSize',5);
box off;

xlabel('Norm. shortening velocity (v/vmax)');
ylabel('Norm. force (f/fiso)');





