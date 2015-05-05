%% 分割点位置与最坏形状的关系 (a为待分割边,abc在上,ade在下)
a=1;
b=0.92;
c=0.52;
d=0.82;
e=0.32;
% divide a
a1=0.01:0.01:0.99;
a2=a-a1;
f=sqrt(a1.^2+b^2-2*a1*b*(a^2+b^2-c^2)/(2*a*b));
g=sqrt(a1.^2+d^2-2*a1*d*(a^2+d^2-e^2)/(2*a*d));
worstShape=min(min(min(shapeCoeff(b,f,a1),shapeCoeff(d,g,a1)),shapeCoeff(c,f,a2)),shapeCoeff(e,g,a2));
plot(a1,worstShape);

% 事实证明是单峰函数，可以用三分法求极大值