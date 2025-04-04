alpha = linspace(0, 2*pi, 100);
beta = linspace(0, pi, 100);
[Alpha, Beta] = meshgrid(alpha, beta);
x = sin(Beta).*cos(Alpha);
y = sin(Beta).*sin(Alpha);
z = cos(Beta);
figure;
surf(x, y, z, 'FaceAlpha', 0.5);
axis equal;
title('3D Convex Example');
xlabel('x(1)');
ylabel('x(2)');
zlabel('x(3)');
hold on
f1 = @(x) x(1);
f2 = @(x) x(2);
f3 = @(x) x(3);
x0 = [0;0;0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
c = 'x(1).^2 + x(2).^2 + x(3).^2 -1';
ceq = '[]';
nonlcon = str2func(['@(x) deal(' c ',' ceq ')']);
[x]= MOP_3_obj_ws(f1, f2, f3, x0, A, b, Aeq, beq, lb, ub, nonlcon);
function [x] = MOP_3_obj_ws(f1, f2, f3, x0, A, b, Aeq, beq, lb, ub, nonlcon)
    eps1 = 0;
    while eps1 < 1 + 0.01
        eps2 = 0;
        while eps2 < 1-eps1 +0.01
            f = @(x) eps1*f1(x) + eps2*f2(x) + (1 - eps1 - eps2)*f3(x);
            [x]=fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon);
            plot3(f1(x), f2(x), f3(x), '*');
            hold on
            eps2 = eps2 + 0.01;
        end
        eps1 = eps1 + 0.01;
    end
end