A= randi(3,5,4)
rank(A)
if rank(A) ~= 4
    error('Nem teljes oszlopszar')
end
c = randi(12,5,1)
rref([A,c])
pinv(A)*c
B=transpose(A)
b=randi(8,4,1)
pinv(B)*b
eye(5)-pinv(B)*B
