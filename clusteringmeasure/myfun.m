function f=myfun(G_a,X,W,mu,V)
syms G
G_b=G*mu;
for i=1:V
   G_b=G_b+(X{i}')*(W{i}.*(X{i}*G));
end
f=matlabFunction(G_a+G_b);
end
