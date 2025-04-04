x_star = cell(2, 1);
y_ideal = zeros(2,1);
f_function = cell(1, 2);
f_string = cell(1, 2); 
for i = 1:2
    f_string{i} = input(['Enter objective function number ' num2str(i) ': '], 's');
    f_function{i} = str2func(['@(x) ' f_string{i}]);
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
x0=input('Enter phase x0: ');
new_ceq = cell(1, 2);
new_nonlcon1 = cell(1, 2);
for i = 1:2
    [~, y_ideal(i)]=fmincon(f_function{i},x0,A,b,Aeq,beq,lb,ub,nonlcon);
    new_ceq{i} = [f_string{i} '-' num2str(y_ideal(i))];
    new_nonlcon1{3-i} = str2func(['@(x) deal(' c ',[' ceq ';' new_ceq{i} '])']);
end
y_nadir = zeros(2,1);
for i = 1:2
    nonlcon = new_nonlcon1{i};
    [~, y_nadir(i)]=fmincon(f_function{i},x0,A,b,Aeq,beq,lb,ub,nonlcon);
end
disp("Ideal point: ");
disp(y_ideal);
disp("Nadir point: ");
disp(y_nadir);