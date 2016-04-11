%% Matematika 7. h�zi feladat
% Zink Martin OZCBJ8
% (Octave, Matlab vagy Sage haszn�lat�val) Vegyen egy (saj�t) pixelform�tum� sz�rke�rnyalatos f�nyk�pet �s k�sz�tsen bel�le egy m�trixot, melynek minden eleme a k�p egy pontj�nak �rnyalalt�t adja meg (pl. Octave-ban �s Matlabban az imread/imwrite paranccsal lehet beolvasni/ki�rni, de van parancs a sz�nesb�l sz�rk�re val� konverzi�ra is). Legyen e m�trix rangja r. �rja fel a m�trix szingul�ris felbont�s�t, majd abb�l k�sz�tsen k�t olyan k�p-m�trixot, melyek egyik�nek 10, m�sik�nak r/2 (eg�sz r�sze) a rangja, �s ,,legk�zelebb'' van az eredeti m�trixhoz (ez �gy sz�molhat�, hogy az els� 10, illetve els� r/2 tagj�t adja �ssze a szingul�ris felbont�s diadikus alakj�nak). Jelen�tse meg e k�t k�pet az eredeti mellett. Ezut�n mindk�t k�pre sz�m�tsa ki az elhagyott szingul�ris �rt�kek n�gyzet�sszeg�nek gy�k�t (ez fogja m�rni a k�t m�trix t�vols�g�t Frobenius-norm�ban), �s adja meg a legnagyobb elhagyott szingul�ris �rt�ket (ez fogja m�rni a k�t m�trix t�vols�g�t 2-norm�ban � a norma fogalma a k�vetkez� el�ad�s anyaga, akkor fogjuk tiszt�zni ezek matematikai tartalm�t). (Akinek van kedve, �s sz�nes nyomtat�ja, megteheti hogy a h�rom sz�n�sszetev� mindegyik�re elv�gzi a szingul�ris felbont�st, majd a fenti elj�r�st, �gy h�rom sz�nes k�pet jelen�t meg).	(2 pont) 
% Hat�rid�: 2015-04-12 10:15
clear all
close all
imgname = 'billy_grey'; %original picture 
type = '.bmp';  %picture ext(i.e .png, .bmp)

%% K�p beolvas�sa

A = double(imread(strcat(imgname,type)));
[U,S,V] = svd(A);
r = rank(A)
S10 = [[S(1:10,1:10),zeros(10,r-10)];zeros(r-10,r)];
SR2 = [[S(1:r/2,1:r/2),zeros(r/2,r-r/2)];zeros(r-r/2,r)];
B = U*S10*V';
C = U*SR2*V';

%%
%K�p ki�r�sa
imwrite(uint8(B),strcat(imgname,strcat('_10',type)));
imwrite(uint8(C),strcat(imgname,strcat('_R2',type)));
%%
%BFdistance: az els� k�p n�gyzet�sszeg�nek gy�ke (ez fogja m�rni a k�t
%m�trix t�vols�g�t Frobenius-norm�ban)
%
%BFmax az els� k�p legnagyobb elhagyott szingul�ris �rt�ke (ez fogja m�rni
%a k�t m�trix t�vols�g�t 2-norm�ban)
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
%CFdistance: a m�sodikk�p n�gyzet�sszeg�nek gy�ke (ez fogja m�rni a k�t
%m�trix t�vols�g�t Frobenius-norm�ban)
%
%CFmax: a m�sodik k�p legnagyobb elhagyott szingul�ris �rt�ke (ez fogja
%m�rni a k�t m�trix t�vols�g�t 2-norm�ban)
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