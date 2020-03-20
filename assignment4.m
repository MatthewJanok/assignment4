%Part 3

clear
clc

% NumElec = 200;
nx = 20;
ny = 40;
Vx = 0.8;
Vy = 0;
delx = 0.05e-7;
dely = 0.05e-7;
Acond = 1;
Bcond = 10E-2;
oldCurr = 0;

for m = 1:4

G = sparse(nx*ny,nx*ny);
Vv = zeros(nx*ny,1);
V = zeros(nx,ny);
B = zeros((nx*ny),1);

x = linspace(0,2e-7,40);
y = linspace(0,1e-7,20);
[X,Y] = meshgrid(x,y);


%Part 2 from assignment 2

cMap = zeros(nx, ny);

for j = 1:ny
    for i = 1:nx
        cMap(i,j) = Acond;
        if ((j<((2/3)*ny))&&((j>((1/3)*ny))&&(i<((m/9)*nx))))|| ((j<((2/3)*ny))&&((j>((1/3)*ny))&&(i>(nx-(m/9)*nx))))
            cMap(i,j) = Bcond;
        end
    end
end
figure(17)
surf(X,Y,cMap)    
title('Conduction Map')
xlabel('ny Value')
ylabel('nx Value')

%Initialize Left Boundary Conditions
for i = 1:nx*ny
    B(i,1) = 0;
    B(i,ny) = 0;
end
        
%Set diagonal
for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        %Set the Boundary Nodes
        if j == 1
            G(n,:) = 0;
            G(n,n) = 1; 
%             B(n,1) = 0;
            
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = Vx;
            
        elseif i == 1
            %Mapping
            nxm = (i)+(j-2)*nx;
            nxp = (i)+(j)*nx;
            nyp = (i+1)+(j-1)*nx;
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 0;
            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);            
            ryp = ((cMap(i,j) + cMap(i+1,j))/2);

            G(n,n) = -(rxp+rxm+ryp);
            G(n,nxm) = rxm;
            G(n,nyp) = ryp;
            G(n,nxp) = rxp;
            
        elseif i == nx
            %Mapping
            nym = (i-1)+(j-1)*nx;
            nxp = (i)+(j)*nx;
            nxm = (i)+(j-2)*nx;
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 0;
            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);            
            rym = ((cMap(i,j) + cMap(i-1,j))/2);

            G(n,n) = -(rxm+rym+rxp);
            G(n,nym) = rym;
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;            

        else
            %Mapping
            nym = (i-1)+(j-1)*nx;
            nyp = (i+1)+(j-1)*nx;
            nxm = (i)+(j-2)*nx;
            nxp = (i)+(j)*nx;

            rym = ((cMap(i,j) + cMap(i-1,j))/2);
            ryp = ((cMap(i,j) + cMap(i+1,j))/2);            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);        
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);       
      
   
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            G(n,nxm) = rxm;
            G(n,nxp) = rxp; 
            B(n,1) = 0;
            


        end
   
    end
end

Vv = G\B;

for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        V(i,j) = Vv(n,1);
    end
end

figure(18)
surf(X,Y,V)
title('V Plot')
xlabel('ny Value')
ylabel('nx Value')



%Gradient of V
for j = 1:ny
    for i = 1:nx
        if j == 1
            Ex(i,j) = (V(i,j+1)-V(i,j));
        elseif j == ny
            Ex(i,j) = (V(i,j)-V(i,j-1));
        else
            Ex(i,j) = (V(i,j+1)-V(i,j-1))*0.5;
        end
        if i == 1
            Ey(i,j) = (V(i+1,j)-V(i,j));
        elseif i == nx
            Ey(i,j) = (V(i,j)-V(i-1,j));
        else
            Ey(i,j) = (V(i+1,j)-V(i-1,j))*0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;
% figure(6)
% surf(X,Y,Ex)
% title('Ex Field')
% xlabel('ny Value')
% ylabel('nx Value')

% figure(7)
% surf(X,Y,Ey)
% title('Ey Field')
% xlabel('ny Value')
% ylabel('nx Value')

E = sqrt(Ex.^2+Ey.^2);
figure(19)
surf(X,Y,E)
title('Total E Field')
xlabel('ny Value')
ylabel('nx Value')

eflowx = cMap.*-Ex;
eflowy = cMap.*-Ey;

% figure(9)
% surf(X,Y,eflowx)
% title('Eflow X')
% xlabel('ny Value')
% ylabel('nx Value')
% 
% figure(10)
% surf(X,Y,eflowy)
% title('Eflow Y')
% xlabel('ny Value')
% ylabel('nx Value')


eflow = sqrt(eflowx.^2 + eflowy.^2);

% figure(11)
% surf(X,Y,eflow)
% title('Eflow Total')
% xlabel('ny Value')
% ylabel('nx Value')

C0 = sum(eflowx(1,:));
Cnx = sum(eflowx(nx,:));
Curr = (C0+Cnx)*0.5;

figure(20)
quiver(X,Y,Ex,Ey)
axis ([0 200e-9 0 100e-9])

figure(21)
hold on
plot([(m-1)*2 (m)*2],[oldCurr Curr], 'k')
title('Average Current at Different Bottleneck Widths')
xlabel('Bottleneck Portion of 10')
ylabel('Current (A)')
hold off

oldCurr = Curr;

pause(1)
end
