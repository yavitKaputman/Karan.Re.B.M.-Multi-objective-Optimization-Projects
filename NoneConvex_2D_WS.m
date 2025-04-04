f1 = @(x) -3*x(1)-2*x(2)+3;
f2 = @(x) -x(1)-3*x(2)+1;
x0 = [0.5;0.5];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0, 0];
ub = [];
c = '(x(1)-1).^3 + x(2)';
ceq = '[]';
nonlcon = str2func(['@(x) deal(' c ',' ceq ')']);
[x]= MOP_2_obj_ws(f1, f2, x0, A, b, Aeq, beq, lb, ub, nonlcon);
function [x] = MOP_2_obj_ws(f1, f2, x0, A, b, Aeq, beq, lb, ub, nonlcon)
    eps = 0;
    while eps < 1 + 0.01
        f = @(x) eps*f1(x) + (1-eps)*f2(x);
        [x]=fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon);
        plot(f1(x), f2(x), '*')
        hold on
        eps = eps + 0.01;
    end
    xlim([0, 2]);
    ylim([-2, 1]);
    title('2D None-Convex Example');
    xlabel('x(1)');
    ylabel('x(2)');
    hold on
end