%% Méréselmélet 1. házifeladat
% Zink Martin OZCBJ8
% 2015/2016/2 félév

%%  1.1. feladat
% Állítson elõ u(n), n=0, 1, … diszkrét értéksorozatot multi-szinusz
% generátor „segítségével” mégpedig úgy, hogy diszkrét jel M harmonikus
% komponensbõl álljon, az egyes harmonikus komponensek amplitúdója
% egységnyi, kezdõfázisa véletlen, és a sorozat várható értéke pedig 
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
% Készítsen multi-szinusz analizátort, amely az u(n) sorozatból kiszámítja az egyes harmonikus komponensek
% amplitúdóját és fázisát! Ellenõrzésképpen végezze is el a számítást (max. 2 pont)!
FFU = fft(u,N);             % FFT
FFA=2*abs(FFU(1:N/2+1))/N;  % Amplitude
inputa=FFA(1:M);
%%
% Harmonikusok amplitúdói
stem(FFA,'fill')
xlim([0 500])
title('Single-Sided Amplitude Spectrum')
xlabel('f')
ylabel('|U(f)|')


%%
% Harmonikusok fázisai
figure
FFP=angle(FFU(1:N/2+1)); % Phase
inputp=FFP(1:M);
stem(FFA.*FFP,'fill')
xlim([0 500])
title('Single-Sided Phase Spectrum')
xlabel('f')
%%
% Az amplitúdó szemmel láthatóan helyes, míg a fáziskarakterika
% összehasonlítható az eredeti phase vector adataival

%% 1.3. feladat
% Ez az u(n) sorozat legyen a bemenõjele a modellezendõ/adaptálandó
% rendszernek, melynek átviteli függvénye
%
% $$ B = \frac{(1-r)z^{-1}}{1-rz^{-2}} $$
%
% A multi-szinusz analizátorral mérje meg a modellezendõ/adaptálandó rendszer kimenõjelét, és a multi-szinusz
% frekvenciákon határozza meg az átvitel abszolút értékét és fázisát (max. 2 pont)!

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


FFP=angle(FFY(1:N/2+1)); % fázis
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
% Amplitúdó karakterisztika
stem(ha(2:end),'fill')
title('Amplitudo atvitel')
xlabel('f')
ylabel('|H(f)|')
%%
% Fázis karakterisztika
figure
stem(hp,'fill')
title('Fazis atvitel')
xlabel('f')
ylabel('p(f)')

%% 1.4. feladat
% Az adaptálandó rendszert lineáris kombinátorral igyekszünk modellezni. Ennek kimenõjele
% TODO equation
% Juttassa el a modellezendõ/adaptálandó rendszert az állandósult állapotáig, majd alkalmas eljárással határozza
% meg (2) súlyozó együtthatóit! Mutassa be, hogy milyen eljárást választott! (A lineáris kombinátor súlytényezõibõl
% alkotott állóvektort W, a regressziós vektort X jelölje!) (max. 3 pont)


N=20; % együtthatók száma (P)
% Vizsgálandó mintasor kivágása. Állandósult állapot elérése után M minta.
offset = 5000;
u_win = u(offset:offset+M-1);
y_win = y(offset:offset+M-1);
% A lineáris kombinátor felfogható egy P hosszu FIR szüröként, aminek a bemeneti vektora X.
X = zeros(N,M);
for i=1:N
    X(i,i:end) = u_win(1:M+1-i);
end
% X autokorrelációs mátrixának (R) számítása. A várható érték képzés átlagolásra egyszerüsödik.
R = zeros(N);
for i=1:M
    R = R + X(:,i)*X(:,i)';
end
R = R/M;
% Hasonló módon X és y keresztkorrelációjának (P) számítása.
P = zeros(N,1);
for i=1:M
    P = P + X(:,i)*y_win(i);
end
P = P/M;
% Lineáris kombinátor optimális együtthatóinak (W) számolása
W = inv(R)*P

y_win_prime = filter(W,1,u_win); % becsült kimenet
e = y_win - y_win_prime; % hiba számítása

%% 1.5. feladat
% Fejtse sorba az (1) átviteli függvényt! Vesse össze a sorfejtett alak és a lineáris kombinátor együtthatóit az (1)
% átviteli függvényû rendszer súlyfüggvényével! Ábrázolja és vesse össze a modellezendõ rendszer és a lineáris
% kombinátor - méréssel meghatározott - átvitelének abszolút értékét és fázisát lineáris skálákat alkalmazva
% (max. 3 pont)!

% ?????

%% 2.A. feladat
% Próbálja ki az ?-LMS eljárást az 1. feladat esetére:
%
%
% A paraméterek nulla kezdeti értékébõl indulva futtassa az algoritmust a (közelítõ) megoldás megtalálásáig. Ezt
% követõen r értékét csökkentse q-val, majd folytassa a futtatást az új megoldás megtalálásáig. A bátorsági
% tényezõt Ön válassza meg! Indokolja választását (max. 2 pont)! Rajzolja ki az együtthatók alakulását az iterációs
% lépések függvényében (konvergencia diagram). Beadandó a program kommentezett listája (max. 2 pont) és az
% együtthatók konvergencia diagramja (max. 4 pont)!

q=0.21;
N = 20;
LENGTH=2^14;
% A vizsgálandó módosított rendszer
num_m = [0 (1-(r-q)) 0];
den_m = [1 0 -(r-q)];
% Módosított rendszer kimenetének elöállítása a gerjesztés hatására
y_m = filter(num_m,den_m,u);
% A két rendszer válaszának összerakása.
y_LMS = [y(1:floor(LENGTH/2)) y_m(floor(LENGTH/2)+1:end)];
% Lineáris kombinátor gerjesztésének elöállítása teljes mintaszámra
X = zeros(N,LENGTH);
for i=1:N
    X(i,i:end) = u(1:LENGTH+1-i);
end

alpha = 1/(N*u*u'/LENGTH)

%%
% Bátorsági tényezõ ( $$ \mu $$  )meghatározása.
% Minél kisebb $$ \mu $$, annál pontosabb lesz a megtalált közelítés,
% viszont annál lassabb lesz az algoritmus.
% Éppen ezért minden alkalmazáshoz megkell határozni az elfogadható
% kompromisszumot.

lambda = max(eig(R));
mu = 1/(2*lambda)/200
%%
W = zeros(N,1);
W_backup = zeros(N,LENGTH);
% Lépésenként közelítjük
e = zeros(1,LENGTH);
for i=1:LENGTH
    e(i) = y_LMS(i) - X(:,i)'*W; % hiba számolása
    W = W + (alpha*e(i))/(X(:,i)'*X(:,i))*X(:,i); % a hiba információjából új W számolása
    W_backup(:,i) = W; % új W eltárolása
end

figure(10)
hold all
for i=1:N
    plot(W_backup(i,:))
end
title('Konvergencia diagram LMS alpha algoritmushoz');

%% 3. feladat
% Végezze el a modellillesztést mindkét átviteli függvényre (r és r-q esetek) az
% alakú, végtelen impulzusválaszú modell alkalmazásával is! Az illesztés során a véges impulzusválaszú
% problémára visszavezetés módszerét (equation-error formulation) alkalmazza! Beadandó a program
% kommentezett listája (max. 2 pont) és az együtthatók konvergencia diagramja (max. 4 pont)!

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