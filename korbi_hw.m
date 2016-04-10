%% M�r�selm�let 1. h�zifeladat
% Zink Martin OZCBJ8
% 2015/2016/2 f�l�v

%%  1.1. feladat
% �ll�tson el� u(n), n=0, 1, � diszkr�t �rt�ksorozatot multi-szinusz
% gener�tor �seg�ts�g�vel� m�gpedig �gy, hogy diszkr�t jel M harmonikus
% komponensb�l �lljon, az egyes harmonikus komponensek amplit�d�ja
% egys�gnyi, kezd�f�zisa v�letlen, �s a sorozat v�rhat� �rt�ke pedig 
% nulla legyen (max. 2 pont)!


close all
clear all

N = 2^14;               % No Samples
M = 230;                % No Harmonics

phase = 2*pi*rand(M,1)-pi;     % Random phase

u = zeros([1 N]);

for i=1:M            % Adding harmonics
    u = u+cos(i*2*pi*(0:1:N-1)/N+phase(i));
end

%% 1.2. feladat
% K�sz�tsen multi-szinusz analiz�tort, amely az u(n) sorozatb�l kisz�m�tja az egyes harmonikus komponensek
% amplit�d�j�t �s f�zis�t! Ellen�rz�sk�ppen v�gezze is el a sz�m�t�st (max. 2 pont)!
FFU = fft(u,N);             % FFT
FFA=2*abs(FFU(1:N/2+1))/N;  % Amplitude
inputa=FFA(1:M);
%%
% Harmonikusok amplit�d�i
stem(FFA,'fill')
xlim([0 500])
title('Single-Sided Amplitude Spectrum')
xlabel('f')
ylabel('|U(f)|')


%%
% Harmonikusok f�zisai
figure
FFP=angle(FFU(1:N/2+1)); % Phase
inputp=FFP(1:M);
stem(FFA.*FFP,'fill')
xlim([0 500])
title('Single-Sided Phase Spectrum')
xlabel('f')
%%
% Az amplit�d� szemmel l�that�an helyes, m�g a f�ziskarakterika
% �sszehasonl�that� az eredeti phase vector adataival

%% 1.3. feladat
% Ez az u(n) sorozat legyen a bemen�jele a modellezend�/adapt�land�
% rendszernek, melynek �tviteli f�ggv�nye
%
% $$ B = \frac{(1-r)z^{-1}}{1-rz^{-2}} $$
%
% A multi-szinusz analiz�torral m�rje meg a modellezend�/adapt�land� rendszer kimen�jel�t, �s a multi-szinusz
% frekvenci�kon hat�rozza meg az �tvitel abszol�t �rt�k�t �s f�zis�t (max. 2 pont)!

figure
r=0.80;
% System
num = [0 (1-r) 0];
den = [1 0 -r];
% System output from s input
y = filter(num,den,u);


FFY = fft(y,N);    % FFT
FFA=2*abs(FFY(1:N/2+1))/N; % Amplitude
outputa=FFA(1:M);
stem(FFA,'fill')
xlim([0 500])
title('Single-Sided Amplitude Spectrum')
xlabel('f')
ylabel('|Y(f)|')


FFP=angle(FFY(1:N/2+1)); % f�zis
outputp=FFP(1:M);
% stem(FFP,'fill')
% xlim([0 500])
% title('Single-Sided Phase Spectrum')
% xlabel('f')

for i=1:M
    ha(i)=outputa(i)/inputa(i); % Abs
    hp(i)=outputp(i)-inputp(i); % Arg
end

%%
% Amplit�d� karakterisztika
stem(ha(2:end),'fill')
title('Amplitudo atvitel')
xlabel('f')
ylabel('|H(f)|')
%%
% F�zis karakterisztika
figure
stem(hp,'fill')
title('Fazis atvitel')
xlabel('f')
ylabel('p(f)')

%% 1.4. feladat
% Az adapt�land� rendszert line�ris kombin�torral igyeksz�nk modellezni. Ennek kimen�jele
% TODO equation
% Juttassa el a modellezend�/adapt�land� rendszert az �lland�sult �llapot�ig, majd alkalmas elj�r�ssal hat�rozza
% meg (2) s�lyoz� egy�tthat�it! Mutassa be, hogy milyen elj�r�st v�lasztott! (A line�ris kombin�tor s�lyt�nyez�ib�l
% alkotott �ll�vektort W, a regresszi�s vektort X jel�lje!) (max. 3 pont)


N=20; % egy�tthat�k sz�ma (P)
% Vizsg�land� mintasor kiv�g�sa. �lland�sult �llapot el�r�se ut�n M minta.
offset = 5000;
u_win = u(offset:offset+M-1);
y_win = y(offset:offset+M-1);
% A line�ris kombin�tor felfoghat� egy P hosszu FIR sz�r�k�nt, aminek a bemeneti vektora X.
X = zeros(N,M);
for i=1:N
    X(i,i:end) = u_win(1:M+1-i);
end
% X autokorrel�ci�s m�trix�nak (R) sz�m�t�sa. A v�rhat� �rt�k k�pz�s �tlagol�sra egyszer�s�dik.
R = zeros(N);
for i=1:M
    R = R + X(:,i)*X(:,i)';
end
R = R/M;
% Hasonl� m�don X �s y keresztkorrel�ci�j�nak (P) sz�m�t�sa.
P = zeros(N,1);
for i=1:M
    P = P + X(:,i)*y_win(i);
end
P = P/M;
% Line�ris kombin�tor optim�lis egy�tthat�inak (W) sz�mol�sa
W = inv(R)*P

y_win_prime = filter(W,1,u_win); % becs�lt kimenet
e = y_win - y_win_prime; % hiba sz�m�t�sa

%% 1.5. feladat
% Fejtse sorba az (1) �tviteli f�ggv�nyt! Vesse �ssze a sorfejtett alak �s a line�ris kombin�tor egy�tthat�it az (1)
% �tviteli f�ggv�ny� rendszer s�lyf�ggv�ny�vel! �br�zolja �s vesse �ssze a modellezend� rendszer �s a line�ris
% kombin�tor - m�r�ssel meghat�rozott - �tvitel�nek abszol�t �rt�k�t �s f�zis�t line�ris sk�l�kat alkalmazva
% (max. 3 pont)!

% ?????

%% 2.A. feladat
% Pr�b�lja ki az ?-LMS elj�r�st az 1. feladat eset�re:
%
%
% A param�terek nulla kezdeti �rt�k�b�l indulva futtassa az algoritmust a (k�zel�t�) megold�s megtal�l�s�ig. Ezt
% k�vet�en r �rt�k�t cs�kkentse q-val, majd folytassa a futtat�st az �j megold�s megtal�l�s�ig. A b�tors�gi
% t�nyez�t �n v�lassza meg! Indokolja v�laszt�s�t (max. 2 pont)! Rajzolja ki az egy�tthat�k alakul�s�t az iter�ci�s
% l�p�sek f�ggv�ny�ben (konvergencia diagram). Beadand� a program kommentezett list�ja (max. 2 pont) �s az
% egy�tthat�k konvergencia diagramja (max. 4 pont)!

q=0.21;
N = 20;
LENGTH=2^14;
% A vizsg�land� m�dos�tott rendszer
num_m = [0 (1-(r-q)) 0];
den_m = [1 0 -(r-q)];
% M�dos�tott rendszer kimenet�nek el��ll�t�sa a gerjeszt�s hat�s�ra
y_m = filter(num_m,den_m,u);
% A k�t rendszer v�lasz�nak �sszerak�sa.
y_LMS = [y(1:floor(LENGTH/2)) y_m(floor(LENGTH/2)+1:end)];
% Line�ris kombin�tor gerjeszt�s�nek el��ll�t�sa teljes mintasz�mra
X = zeros(N,LENGTH);
for i=1:N
    X(i,i:end) = u(1:LENGTH+1-i);
end

alpha = 1/(N*u*u'/LENGTH)

%%
% B�tors�gi t�nyez� ( $$ \mu $$  )meghat�roz�sa.
% Min�l kisebb $$ \mu $$, ann�l pontosabb lesz a megtal�lt k�zel�t�s,
% viszont ann�l lassabb lesz az algoritmus.
% �ppen ez�rt minden alkalmaz�shoz megkell hat�rozni az elfogadhat�
% kompromisszumot.

lambda = max(eig(R));
mu = 1/(2*lambda)/200
%%
W = zeros(N,1);
W_backup = zeros(N,LENGTH);
% L�p�senk�nt k�zel�tj�k
e = zeros(1,LENGTH);
for i=1:LENGTH
    e(i) = y_LMS(i) - X(:,i)'*W; % hiba sz�mol�sa
    W = W + (alpha*e(i))/(X(:,i)'*X(:,i))*X(:,i); % a hiba inform�ci�j�b�l �j W sz�mol�sa
    W_backup(:,i) = W; % �j W elt�rol�sa
end

figure(10)
hold all
for i=1:N
    plot(W_backup(i,:))
end
title('Konvergencia diagram LMS alpha algoritmushoz');

%% 3. feladat
% V�gezze el a modellilleszt�st mindk�t �tviteli f�ggv�nyre (r �s r-q esetek) az
% alak�, v�gtelen impulzusv�lasz� modell alkalmaz�s�val is! Az illeszt�s sor�n a v�ges impulzusv�lasz�
% probl�m�ra visszavezet�s m�dszer�t (equation-error formulation) alkalmazza! Beadand� a program
% kommentezett list�ja (max. 2 pont) �s az egy�tthat�k konvergencia diagramja (max. 4 pont)!

Ri = inv(R);
W2 = zeros(N,1);
W2_backup = zeros(N,LENGTH);
e2 = zeros (1,LENGTH);
for i = 1:LENGTH
    e2(i) = y_LMS(i) - X(:,i)'*W;
    W2 = W2 + 2*mu*Ri(:,i)*X(:,i)'*e2(i)
    W2_backup(:,i) = W;
end

figure(11)
hold all
for i=1:N
    plot(W_backup(i,:))
end
title('Konvergencia diagram LMS alpha algoritmushoz');