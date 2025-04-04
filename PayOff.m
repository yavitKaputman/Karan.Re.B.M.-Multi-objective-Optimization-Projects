p = input('Enter the number of objective functions: ');
f_function = cell(1, p);
pay_off = zeros(p, p);
y_ideal = zeros(p, 1);
for i = 1:p
    f_i = input(['Enter objective function number ' num2str(i) ': '], 's');
    f_function{i} = str2func(['@(x) ' f_i]);
end 
A = input('Enter A: ');
b = input('Enter b: ');
Aeq = input('Enter Aeq: ');
beq = input('Enter beq: ');
lb = input('Enter lb: ');
ub = input('Enter ub: ');
c = input('Enter c: ', 's');
ceq = input('Enter ceq: ', 's');
nonlcon = str2func(['@(x) deal(' c ',' ceq ')']);
x0=input('Enter x0: ');
for i = 1:p
    [x, y_ideal(i)]=fmincon(f_function{i},x0,A,b,Aeq,beq,lb,ub,nonlcon);
    for j = 1:p
        pay_off(j, i) = f_function{j}(x);
    end
end
y_nadir = max(pay_off, [], 2);
disp("Ideal point: ");
disp(y_ideal);
disp("pay-off matrix: ");
disp(pay_off);
disp("Estimated Nadir point according to pay-off method: ");
disp(y_nadir);