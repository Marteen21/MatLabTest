% W = [1;1;0]
% V = [0 1 0;2 0 0]
% VIO = [0 1 0;2 0 0;1 2 0]
% VIW = [V W]
% P = VIO*inv(VIW)

%%%% DansGame
clear
v1=randi(3,3,1)
v2=randi(3,3,1)
v3=randi(3,3,1)

V = [v1,v2]
W = [v3]
U = [V,W]
if rank(U)~=3
    error('restart m8, the rank is not quite right sry m8')
else
    P= [V,zeros(3,1)]*inv(U)
end
%%% Checking
E = P*randi(3,3,1)

rank([v1,v2,E])