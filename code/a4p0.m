% Jay Maini 101037537
clear all
close all
set(0, 'DefaultFigureWindowStyle', 'docked')

%po = linspace(0.05, 0.95, 19);
v0 = linspace(0.1,10,20);
figure
for i=1:length(v0)
    Cur(i) = CurrentFinder(v0(i),0.5);
end
plot(Cur, v0)
hold on
ylabel('Input Voltage')
xlabel('Average Current Density (A/nm^2)')
title('Current versus Voltage Characteristic')
P = polyfit(Cur, v0, 1)


function J_avg = CurrentFinder(v0, percent_open)
    set(0, 'DefaultFigureWindowStyle','docked')

    %v0 = 1;
    %unequal array size
    nx = 200;
    ny = 100;
    V = zeros(nx,ny);
    G = sparse(nx*ny, nx*ny);
    F = zeros(nx*ny,1);

    %Easy adjustment for bottleneck size - useful for function
    %percent_open = 0.33;
    percent_boundary = (1 - percent_open)/2;
    b_height_1 = round(percent_boundary*ny);
    b_height_2 = round(ny-percent_boundary*ny);

    %Make conductivity array 
    cMap = zeros(nx,ny) + 1;
    %Create bottleneck areas
    BNx = [(2/5)*nx (3/5)*nx];
    %BNx = [b_height_1 b_height_2];
    BNy1 = [1 b_height_1];
    BNy2 = [b_height_2 ny];
    %Modify cMap values

    for o = BNx(1):BNx(2)
        for p = BNy1(1):BNy1(2)
            cMap(o,p) = 1e-2;
        end
        for u = BNy2(1):BNy2(2)
            cMap(o,u) = 1e-2;
        end
    end
    
    for i = 1:nx
        for j = 1:ny
            m = j+(i-1)*ny;
            %Edge cases (BC's nodes)
            if (i == 1)
                F(m) = v0;
                G(m,m) = 1;
            elseif (i == nx)
                F(m) = 0;
                G(m,m) = 1;
            elseif (j == 1)
                %Bottom side - only three resistors
                npy = j+1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;
                nmy = j-1 + (i-1)*ny;

                %Calculate resistances based on conductivity map
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rpy = (cMap(i,j) + cMap(i,j+1))/2;

                G(m,m) = -(rpx+rmx+rpy);
                G(m,npx) = rpx; 
                G(m,nmx) = rmx;
                G(m,npy) = rpy;
            elseif (j == ny)
                %top side - only three resistors
                nmy = j-1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;

                %Calculate resistances based on conductivity map
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rmy = (cMap(i,j) + cMap(i,j-1))/2;

                G(m,m) = -(rpx+rmx+rmy);
                G(m,npx) = rpx; 
                G(m,nmx) = rmx;
                G(m,nmy) = rmy;
            else
                %Inner region - four resistors
                %Define next/previous on the x and y axes
                npy = j+1 + (i-1)*ny;
                nmy = j-1 + (i-1)*ny;
                npx = j + (i-1+1)*ny;
                nmx = j + (i-1-1)*ny;

                %Calculate resistances based on conductivity map
                rpx = (cMap(i,j) + cMap(i+1,j))/2;
                rmx = (cMap(i,j) + cMap(i-1,j))/2;
                rpy = (cMap(i,j) + cMap(i,j+1))/2;
                rmy = (cMap(i,j) + cMap(i,j-1))/2;

                G(m,m) = -(rpx+rmx+rpy+rmy);
                G(m,npx) = rpx; 
                G(m,nmx) = rmx;
                G(m,npy) = rpy;
                G(m,nmy) = rmy;
            end
        end
    end

    %LU Decompsition
    Varr = G\F;
    Varr_temp = reshape(Varr,ny,nx);

    %Output the current obtained
    J = cMap' .*gradient(-Varr_temp);
    %Peak current 
    %J_max = max(J, [], 'all')
    J_avg = mean(J(1,:), 'all');
end




