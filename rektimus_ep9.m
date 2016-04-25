%% Matematika 9. házi feladat
%% (a) feladat
% Generáljunk véletlen elemekkel egy ilyen P mátrixot és egy véletlen
% nemnegatív vv vektort (legyen ?v?1=1?v?1=1)! Milyen pozitivitási
% osztályba tartozik a P mátrix (primitív, reducibilis, irreducibilis)?
clear all;
close all;
n = 12;
for i = 1:n
    pe = rand(1);
    if (i == 1)
        P(i,i+1) = 1-pe;
        P(i,n) = pe;
    elseif (i == 12)
        P(i,1) = 1-pe;
        P(i,n-1) = pe;
    else
        P(i,i-1) = pe;
        P(i,i+1) = 1-pe;
    end
end


sum = 0;
v = rand(1,n)
for i = 1:n
    sum = sum + v(i);
end
v = v/sum;
P
v

%%
% P mátrix irreducibilis, mert a belõle felírt mátrix minden élbõl bejárható

%% (b) feladat
% Vizsgáljuk meg a vPnvPn vektorok viselkedését, ha nn tart végtelenhez!
% Létezik-e e vektorsorozatnak határértéke vagy torlódási pontja?
% Létezik-e a v,vP,…,vPk?1v,vP,…,vPk?1 vektorok átlagának határértéke?
% Mi ezek jelentése?

for i=1:1000    %A végtelen esetet egy igen nagy (1000 db) os esettel közelítettem
    limesN = v*P^i
    if(mod(i,2)==0)
        figure(1)
        title('Páros hatványok');
        plot(limesN)
        hold on
        
    else
        figure(2)
        title('Páratlan hatványok');
        plot(limesN)
        hold on
    end
end

%%
% A vektorsorozat szemmel láthatóan torlódik, határértéke ennek ellenére
% nincsen

%% (c) feladat
r = max(eig(P))
[p_root q] = perron(P,'left')

%------------MATHWORKS PERRON-----------------
[V D]=eig(P)
lambda = diag(D);
rho = max(lambda);
tol = 1e-10; % use some small tolerance
ind = find(lambda>=rho*(1-tol));
V = V(:,ind);
A = -V;
b = zeros(size(V,1),1);
sum = 0;
for i=1:n
    sum = sum + (V(i));
end
Aeq = sum;
beq = 1;
dummy = ones(size(V,2),1);
x = linprog(dummy, A, b, Aeq, beq);
pv = V*x

%---------------------------------------------