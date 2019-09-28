% This script produces an indicator function for an approximated human hand
% in the unit square (0,0)-(1,1)

rotx = @(x,y,c,theta) cos(theta).*(x-c(1)) + sin(theta).*(y-c(2)) + c(1);
roty = @(x,y,c,theta) -sin(theta).*(x-c(1)) + cos(theta).*(y-c(2)) + c(2);

a1=0.32; b1=0.9; p1=0.04; q1=0.23;
th1=0.19*pi; c1=[0.33,0.76];
finger1 = @(x,y) (  ((x-a1)./p1).^2 + ((y-b1)./q1).^2  ) < 1;
finger1r = @(x,y) finger1(rotx(x,y,c1,th1),roty(x,y,c1,th1));
%
a2=0.4; b2=0.9; p2=0.04; q2=0.3;
th2=0.065*pi; c2=[0.41,0.81];
finger2 = @(x,y) (  ((x-a2)./p2).^2 + ((y-b2)./q2).^2  ) < 1;
finger2r = @(x,y) finger2(rotx(x,y,c2,th2),roty(x,y,c2,th2));
%
a3=0.535; b3=0.95; p3=0.05; q3=0.3;
th3=0*pi; c3=[0.3+2/8,0.5];
finger3 = @(x,y) (  ((x-a3)./p3).^2 + ((y-b3)./q3).^2  ) < 1;
finger3r = @(x,y) finger3(rotx(x,y,c3,th3),roty(x,y,c3,th3));
%
a4=0.64; b4=0.9; p4=0.04; q4=0.24;
th4=-0.08*pi; c4=[0.64,0.74];
finger4 = @(x,y) (  ((x-a4)./p4).^2 + ((y-b4)./q4).^2  ) < 1;
finger4r = @(x,y) finger4(rotx(x,y,c4,th4),roty(x,y,c4,th4));
%
a5=0.7; b5=0.5; p5=0.2; q5=0.08;
th5=0.06*pi; c5=[0.71,0.46];
finger5 = @(x,y) (  ((x-a5)./p5).^2 + ((y-b5)./q5).^2  ) < 1;
finger5r = @(x,y) finger5(rotx(x,y,c5,th5),roty(x,y,c5,th5));
%
a6=0.5; b6=0.61; p6=0.25; q6=0.28;
palm = @(x,y) (  ((x-a6)./p6).^2 + ((y-b6)./q6).^2  ) < 1;
%
hand0 = @(x,y) finger1r(x,y) | finger2r(x,y) | finger3r(x,y) | finger4r(x,y)...
              | finger5r(x,y) | palm(x,y);
hand = @(x,y) hand0(x,y);