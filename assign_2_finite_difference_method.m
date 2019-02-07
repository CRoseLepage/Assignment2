% Assignment 2 Finite Difference Method
% Chantel Lepage 100999893

%%Part 1

L=50;
W=50;
Vo=1;


G=sparse(L,W);
V=zeros(1,L);

for i=1:L
    for j=1:W
        n=j+(i-1)*W; % create 4 different conditions one for each side. top and bottom are not set so set to 0
        if i==1 || j==L
            V(n)=0;
            G(:,n)=0;
            G(n,n)=1;
        elseif j==1 || j==ny
            V(n)=0;
            G(:,n)=0;
            G(n,n)=1;
        else
            nxm = j+(i-2)*ny;
            nxp = j+(i)*ny;
            nym = (j-1)+(i-1)*ny;
            nyp = (j+1)+(i-1)*ny;            
            
            G(n,n) = -4;
            G(n,nxm)= 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end        
    end
end





