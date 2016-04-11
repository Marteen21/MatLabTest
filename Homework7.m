%% Matematika 7. házi feladat
% Zink Martin OZCBJ8
% (Octave, Matlab vagy Sage használatával) Vegyen egy (saját) pixelformátumú szürkeárnyalatos fényképet és készítsen belõle egy mátrixot, melynek minden eleme a kép egy pontjának árnyalaltát adja meg (pl. Octave-ban és Matlabban az imread/imwrite paranccsal lehet beolvasni/kiírni, de van parancs a színesbõl szürkére való konverzióra is). Legyen e mátrix rangja r. Írja fel a mátrix szinguláris felbontását, majd abból készítsen két olyan kép-mátrixot, melyek egyikének 10, másikának r/2 (egész része) a rangja, és ,,legközelebb'' van az eredeti mátrixhoz (ez úgy számolható, hogy az elsõ 10, illetve elsõ r/2 tagját adja össze a szinguláris felbontás diadikus alakjának). Jelenítse meg e két képet az eredeti mellett. Ezután mindkét képre számítsa ki az elhagyott szinguláris értékek négyzetösszegének gyökét (ez fogja mérni a két mátrix távolságát Frobenius-normában), és adja meg a legnagyobb elhagyott szinguláris értéket (ez fogja mérni a két mátrix távolságát 2-normában – a norma fogalma a következõ elõadás anyaga, akkor fogjuk tisztázni ezek matematikai tartalmát). (Akinek van kedve, és színes nyomtatója, megteheti hogy a három színösszetevõ mindegyikére elvégzi a szinguláris felbontást, majd a fenti eljárást, így három színes képet jelenít meg).	(2 pont) 
% Határidõ: 2015-04-12 10:15
clear all
close all
imgname = 'billy_grey'; %original picture 
type = '.bmp';  %picture ext(i.e .png, .bmp)

%% Kép beolvasása

A = double(imread(strcat(imgname,type)));
[U,S,V] = svd(A);
r = rank(A)
S10 = [[S(1:10,1:10),zeros(10,r-10)];zeros(r-10,r)];
SR2 = [[S(1:r/2,1:r/2),zeros(r/2,r-r/2)];zeros(r-r/2,r)];
B = U*S10*V';
C = U*SR2*V';

%%
%Kép kiírása
imwrite(uint8(B),strcat(imgname,strcat('_10',type)));
imwrite(uint8(C),strcat(imgname,strcat('_R2',type)));
%%
%BFdistance: az elsõ kép négyzetösszegének gyöke (ez fogja mérni a két
%mátrix távolságát Frobenius-normában)
%
%BFmax az elsõ kép legnagyobb elhagyott szinguláris értéke (ez fogja mérni
%a két mátrix távolságát 2-normában)
BFdistance = 0;
BFmax = S(11,11);
for i = 11:r
    BFdistance =  BFdistance + S(i,i);
    if(S(i,i)>BFmax)
        BFmax = S(i,i);
    end
end
BFdistance = sqrt(BFdistance)
BFmaxs
%%
%CFdistance: a másodikkép négyzetösszegének gyöke (ez fogja mérni a két
%mátrix távolságát Frobenius-normában)
%
%CFmax: a második kép legnagyobb elhagyott szinguláris értéke (ez fogja
%mérni a két mátrix távolságát 2-normában)
CFdistance = 0;
CFmax = S(r/2+1,r/2+1);
for j = r/2+1:r
    CFdistance =  CFdistance + S(j,j);
    if(S(j,j)>CFmax)
        CFmax = S(j,j);
    end
end
CFdistance = sqrt(CFdistance)
CFmax