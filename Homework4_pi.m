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

V = [v1,v2];
W = [v3];
U = [V,W];
X = [v1,v2,v3];
if rank(U)~=3
    error('restart m8, the rank is not quite right')
else
    P = [V,zeros(3,1)]*inv(U)
    Q = [v1,zeros(3,2)]*inv(X)
    %%% Checking
    E = P*randi(3,3,1)
    F = Q*randi(3,3,1)
    if (rank([v1,v2,E])==2)
        disp('Everything seems about right')
    else
        error('Something wrong m8')
    end
        if (rank([v1,F])==1)
        disp('Everything seems about right')
    else
        error('Something wrong m8')
    end
end
