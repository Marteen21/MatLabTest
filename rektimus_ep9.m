%% Matematika 9. h�zi feladat
%% (a) feladat
% Gener�ljunk v�letlen elemekkel egy ilyen P m�trixot �s egy v�letlen 
% nemnegat�v vv vektort (legyen ?v?1=1?v?1=1)! Milyen pozitivit�si
% oszt�lyba tartozik a P m�trix (primit�v, reducibilis, irreducibilis)?
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
sum = 0;
for i = 1:n
    sum = sum + v(i);
end
P
v

%%
% P m�trix irreducibilis, mert a bel�le fel�rt m�trix minden �lb�l bej�rhat�

%% (b) feladat
% Vizsg�ljuk meg a vPnvPn vektorok viselked�s�t, ha nn tart v�gtelenhez! 
% L�tezik-e e vektorsorozatnak hat�r�rt�ke vagy torl�d�si pontja? 
% L�tezik-e a v,vP,�,vPk?1v,vP,�,vPk?1 vektorok �tlag�nak hat�r�rt�ke? 
% Mi ezek jelent�se? 
for i=1:700
    limesN = v*P^i
end
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