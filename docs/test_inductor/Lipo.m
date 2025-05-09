% FEM code (MATLAB) for the example in Chapter 10
% Generated by Wen Ouyang
% The mesh is generated in a straightforward manner for clarity
% Physical model and X,Y notation:
clear all; close all; clc;
% --------------------------------------------------------------(1)
% Parameters
% define BH curve of 1010 steel in Amperes/m.
BH=[ 0 0; 0.2003 238.7; 0.3204 318.3;
    0.40045 358.1; 0.50055 437.7; 0.5606 477.5;
    0.7908 636.6; 0.931 795.8; 1.1014 1114.1;
    1.2016 1273.2; 1.302 1591.5; 1.4028 2228.2;
    1.524 3183.1; 1.626 4774.6; 1.698 6366.2;
    1.73 7957.7; 1.87 15915.5; 1.99 31831;
    2.04 47746.5; 2.07 63662; 2.095 79577.5;
    2.2 159155; 2.4 318310];
% Constants
u0=4*pi*1e-7; % free space permeability
ur_air=1.0000004; % Air relative permeability
ur_cop=0.999991; % Copper relative permeability
acc=0.1; % acceleration factor

Berror=0.01; % error control factor for B
Aerror=0.01; % error control factor for A
Bmax=2; % saturation level limit
inTome=2.54/100; % inch to meter ratio
% Dimension in RELATIVE value to the previous point
X1 = 0.5; % inches, left air region width
X2 = 2; % inches, left tooth width
X3 = 2; % inches, slot width
X4 = 2; % inches, right tooth width
X5 = 0.5; % inches, right air region width
Y1 = 2; % inches, bottom iron height
Y2 = 0.5; % inches, air gap height
Y3 = 2; % inches, coil height
Y4 = 2; % inches, yoke height
% initial values
Binit=0.5; % Tesla, initial value
vinit=interp1(BH(:,1),BH(:,2),Binit,'linear')/Binit
Aboundary=0; % boundary condition set to zero
J=1000; % Amp/in^2, current density
J = J*2^2/(2*2.54/100)^2; %Amp/m^2, current density
status=0; % check if status=1, solution converged
Loop=0; % loop counts
%-------------------------------------------------------------(2)
% Generate meshes and number the nodes of meshes.
dx = 0.5/3; % X direction mesh step, it is a global setting.
% check dx to match it with X1=0.5,X2=X3=X4=2,X5=0.5.
dy1 = 0.4; % Bottom iron Y direction mesh step, note Y1=2.
dy2 = 0.1; % Air gap Y direction mesh step, notice airgap=0.5.
dy3 = 0.4; % Top iron tooth Y direction mesh step, note Y3=2.
dy4 = 1.0; % Top iron yoke Y direction mesh step, note Y4=2.
% Node number information
x1=round(X1/dx+1); x2=round(x1+X2/dx); x3=round(x2+X3/dx);
x4=round(x3+X4/dx); x5=round(x4+X5/dx);
NodeX=x5; % Nodes in X direction
y1=round(Y1/dy1+1); y2=round(y1+Y2/dy2); y3=round(y2+Y3/dy3);
y4=round(y3+Y4/dy4);
NodeY=y4; % Nodes in Y direction
TotalNodes=round(NodeX*NodeY);
TotalMeshes=round((NodeX-1)*(NodeY-1)*2);

% Develop nodal matrix for nodes with actual node coordinates in meters.
% Each node corresponding to each row, node is numbered from bottom to top,
% from left to right.
Nodes=zeros(TotalNodes,2);
% from the bottom to top, from left to right
for i=1:NodeY  % Y-coordinate
    for j=1:NodeX % X-coordinate
        if i<=y1 % Bottom iron region
            Nodes((i-1)*NodeX+j,:)=[(j-1)*dx,(i-1)*dy1]*inTome;
        end
        if (i>y1)&&(i<=y2) % air region
            Nodes((i-1)*NodeX+j,:)=[(j-1)*dx,Y1+(i-y1)*dy2,]*inTome;
        end
        if (i>y2)&&(i<=y3) % Top iron tooth region
            Nodes((i-1)*NodeX+j,:)=[(j-1)*dx,Y1+Y2+(i-y2)*dy3]*inTome;
        end
        if (i>y3) % Top yoke region
            Nodes((i-1)*NodeX+j,:) = [(j-1)*dx; Y1+Y2+Y3+(i-y3)*dy4]*inTome;
        end
    end
end
%--------------------------------------------------------------(3)
% Number the meshes and vertices.
% Assign material number and the initial reluctance.
% Calculate the triangle area.
A=zeros(TotalNodes,1);
B=zeros(TotalMeshes,1);
delta=zeros(TotalMeshes,1);
% Develop mesh element matrix(MEM), each row is for each mesh with
% 3 triangle nodal numbers, material number, and reluctance for
% each element.
% Material 0: Air
% 1: Iron
% 2: Copper
MEM=zeros(TotalMeshes,5);
for i=1:(NodeY-1)
    for j=1:(NodeX-1)
        if i<y1 % air & iron
            if (j<x1)||(j>=x4) % air
                Material=0;
                v=1/(u0*ur_air);
            else % iron
                Material=1;
                v=vinit;
            end
        end

        if (i>=y1)&&(i<y2) % air
            Material=0;
            v=1/(u0*ur_air);
        end
        if (i>=y2)&&(i<y3) % air, iron & copper
            if (j<x1)||(j>=x4) % air
                Material=0;
                v=1/(u0*ur_air);
            elseif (j>=x2)&&(j<x3) % copper
                Material=2;
                v=1/(u0*ur_cop);
            else % iron
                Material=1;
                v=vinit;
            end
        end
        if (i>=y3)&&(i<y4) % air & iron
            if (j<x1)||(j>=x4) % air
                Material=0;
                v=1/(u0*ur_air);
            else % iron
                Material=1;
                v=vinit;
            end
        end
        elemBot=((j-1)*2+1)+((i-1)*2*(NodeX-1));
        % serial number of bottom triangle
        elemTop=elemBot+1; % serial number of top triangle
        vert1=j+(i-1)*NodeX; % Nodal number of vertex at bottom left
        vert2=vert1+1; % Nodal number of vertex at bottom right
        vert3=vert1+NodeX; % Nodal number of vertex at top left
        vert4=vert3+1; % Nodal number of vertex at top right
        MEM(elemBot,:)=[vert1,vert2,vert3,Material,v];
        MEM(elemTop,:)=[vert3,vert2,vert4,Material,v];
        B(elemBot)=Binit; % initial value for element B
        B(elemTop)=Binit;
    end
end
% Calculate the area of triangle, store in delta vector indexed
% by element number
for index=1:TotalMeshes
    i=MEM(index,1); j=MEM(index,2); k=MEM(index,3);
    xi=Nodes(i,1); yi=Nodes(i,2); xj=Nodes(j,1);
    yj=Nodes(j,2); xk=Nodes(k,1); yk=Nodes(k,2);

    area=((xj*yk-xk*yj)+xi*(yj-yk)+yi*(xk-xj))/2;
    delta(index)=area;
end
% --------------------------------------------------------------(4)
% [S] & [T] Matrices and boundary conditions
S=zeros(TotalNodes,TotalNodes);
T=zeros(TotalNodes,1);
for index=1:TotalMeshes
    i=MEM(index,1); j=MEM(index,2); k=MEM(index,3);
    xi=Nodes(i,1); yi=Nodes(i,2); xj=Nodes(j,1);
    yj=Nodes(j,2); xk=Nodes(k,1); yk=Nodes(k,2);
    ai=xj*yk-xk*yj; bi=yj-yk; ci=xk-xj;
    aj=yi*xk-xi*yk; bj=yk-yi; cj=xi-xk;
    ak=xi*yj-xj*yi; bk=yi-yj; ck=xj-xi;
    S(i,i)=S(i,i)+MEM(index,5)*(bi^2+ci^2)/(4*delta(index));
    S(i,j)=S(i,j)+MEM(index,5)*(bi*bj+ci*cj)/(4*delta(index));
    S(i,k)=S(i,k)+MEM(index,5)*(bi*bk+ci*ck)/(4*delta(index));
    S(j,i)=S(i,j);
    S(j,j)=S(j,j)+MEM(index,5)*(bj^2+cj^2)/(4*delta(index));
    S(j,k)=S(j,k)+MEM(index,5)*(bj*bk+cj*ck)/(4*delta(index));
    S(k,i)=S(i,k);
    S(k,j)=S(j,k);
    S(k,k)=S(k,k)+MEM(index,5)*(bk^2+ck^2)/(4*delta(index));
end
% boundary condition, boundaries all set to zero
% bottom and top lines:
for j=1:NodeX
    ZeroNode1=j;
    ZeroNode2=NodeX*(NodeY-1)+j;
    S(ZeroNode1,:)=zeros(1,size(Nodes,1));
    % entire row set to zeros
    S(ZeroNode1,ZeroNode1)=1; % set the diagonal to one
    T(ZeroNode1)=Aboundary; % A=0
    S(ZeroNode2,:)=zeros(1,size(Nodes,1));
    % entire row set to zeros
    S(ZeroNode2,ZeroNode2)=1; % set the diagonal to one
    T(ZeroNode2)=Aboundary; % A=0
end
% right and left lines:
for i=1:NodeY
    ZeroNode1=i*NodeX;

    ZeroNode2=(i-1)*NodeX+1;
    S(ZeroNode1,:)=zeros(1,size(Nodes,1));
    % entire row set to zeros
    S(ZeroNode1,ZeroNode1)=1; % set the diagonal to one
    T(ZeroNode1)=Aboundary; % A=0
    S(ZeroNode2,:)=zeros(1,size(Nodes,1)); % entire row set to zeros
    S(ZeroNode2,ZeroNode2)=1; % set the diagonal to one
    T(ZeroNode2)=Aboundary; % A=0
end
% current density on copper area
for N=1:TotalNodes
    I=0;
    for index=1:TotalMeshes
        i=MEM(index,1);
        j=MEM(index,2);
        k=MEM(index,3);
        if ((N==i)||(N==j)||(N==k))&&(MEM(index,4)==2)
            I=I+J*delta(index);
        end
    end
    if (I~=0)
        T(N)=I/3;
    end
end
%--------------------------------------------------------------(5)
% Calculating until B & A converge
while ~status % continue until convergence
    Loop=Loop+1; Aold=A; Bold=B;
    A=S\T; % calculate new A by inv(S)*T
    %update the reluctivity of mesh elements and calculate the field energy
    for index=1:TotalMeshes
        i=MEM(index,1); j=MEM(index,2); k=MEM(index,3);
        xi=Nodes(i,1); yi=Nodes(i,2); xj=Nodes(j,1);
        yj=Nodes(j,2); xk=Nodes(k,1); yk=Nodes(k,2);
        ai=xj*yk-xk*yj; bi=yj-yk; ci=xk-xj;
        aj=yi*xk-xi*yk; bj=yk-yi; cj=xi-xk;
        ak=xi*yj-xj*yi; bk=yi-yj; ck=xj-xi;
        % update flux density in each element
        Bx=(1/(2*delta(index)))*[ci*A(MEM(index,1))+cj*A(MEM(index,2))+ck*A(MEM(index,3))];
        By=(-1/(2*delta(index)))*[bi*A(MEM(index,1))+bj*A(MEM(index,2))+bk*A(MEM(index,3))];
        B(index)=sqrt(Bx*Bx+By*By);
        if B(index)>=Bmax % saturation control
            B(index)=Bmax;
        end
        % update the reluctivity,v
        if (MEM(index,4)==1) % iron
            vcalc=interp1(BH(:,1),BH(:,2),B(index),'linear')/B(index);
        elseif (MEM(index,4)==2) % copper
            vcalc=1/(ur_cop*u0);
        else
            vcalc=1/(ur_air*u0); % air
        end
%          v=MEM(index,5)+acc*(vcalc-MEM(index,5));
        MEM(index,5)=vcalc;
        % Update [S]
        S=zeros(TotalNodes,TotalNodes);
        for index_1=1:TotalMeshes
            i=MEM(index_1,1); j=MEM(index_1,2); k=MEM(index_1,3);
            xi=Nodes(i,1); yi=Nodes(i,2); xj=Nodes(j,1);
            yj=Nodes(j,2); xk=Nodes(k,1); yk=Nodes(k,2);
            ai=xj*yk-xk*yj; bi=yj-yk; ci=xk-xj;
            aj=yi*xk-xi*yk; bj=yk-yi; cj=xi-xk;
            ak=xi*yj-xj*yi; bk=yi-yj; ck=xj-xi;
            S(i,i)=S(i,i)+MEM(index_1,5)*(bi^2+ci^2)/(4*delta(index_1));
            S(i,j)=S(i,j)+MEM(index_1,5)*(bi*bj+ci*cj)/(4*delta(index_1));
            S(i,k)=S(i,k)+MEM(index_1,5)*(bi*bk+ci*ck)/(4*delta(index_1));
            S(j,i)=S(i,j);
            S(j,j)=S(j,j)+MEM(index_1,5)*(bj^2+cj^2)/(4*delta(index_1));
            S(j,k)=S(j,k)+MEM(index_1,5)*(bj*bk+cj*ck)/(4*delta(index_1));
            S(k,i)=S(i,k);
            S(k,j)=S(j,k);
            S(k,k)=S(k,k)+MEM(index_1,5)*(bk^2+ck^2)/(4*delta(index_1));
        end
        % Energy calculation for each element
        Energy(index)=B(index)^2*delta(index)*MEM(index,5)/2;
    end % end of index loop
    % Reapply boundary condition.
    % left and right lines:
    for i=1:NodeY
        ZeroNode1=i*NodeX;

        ZeroNode2=(i-1)*NodeX+1;
        S(ZeroNode1,:)=zeros(1,size(Nodes,1));
        % entire row set to zeros
        S(ZeroNode1,ZeroNode1)=1; % set the diagonal to one
        T(ZeroNode1)=Aboundary; % A=0
        S(ZeroNode2,:)=zeros(1,size(Nodes,1));
        % entire row set to zeros
        S(ZeroNode2,ZeroNode2)=1; % set the diagonal to one
        T(ZeroNode2)=Aboundary; % A=0
    end
    % bottom and top lines:
    for j=1:NodeX
        ZeroNode1=j;
        ZeroNode2=NodeX*(NodeY-1)+j;
        S(ZeroNode1,:)=zeros(1,size(Nodes,1));
        % entire row set to zeros
        S(ZeroNode1,ZeroNode1)=1; % set the diagonal to one
        T(ZeroNode1)=Aboundary; % A=0
        S(ZeroNode2,:)=zeros(1,size(Nodes,1));
        % entire row set to zeros
        S(ZeroNode2,ZeroNode2)=1; % set the diagonal to one
        T(ZeroNode2)=Aboundary; % A=0
    end
    deltaB=abs(B-Bold); % Check convergence of B
    deltaA=abs(A-Aold); % Check convergence of A
    status=(sqrt(sum(deltaB.^2))/sqrt(sum(B.^2))<=Berror)&(sqrt(sum(deltaA.^2))/sqrt(sum(A.^2))<=Aerror);
    disp('Loop='); disp(Loop);
    if Loop >=50 % Loop control
        break;
    end
end % end of while Loop
%--------------------------------------------------------------(6)
% Draw flux lines, model, meshes
TotalEnergy=sum(Energy);
TotalFlux=max(A);
Az=zeros(NodeY,NodeX); % rearrange the Az distribution in the model
format("default")
Xaxis=zeros(NodeX,1);
Yaxis=zeros(NodeY,1);
for i=1:NodeY

    for j=1:NodeX
        Az(i,j)=A(j+(i-1)*NodeX);
    end
end
% Find x coordinates values for material critical nodes in INCH.
for j=1:NodeX
    if(j<=x1) Xaxis(j)=(j-1)*dx;end
    if(j>x1)&&(j<=x2) Xaxis(j)=(j-x1)*dx+X1;end
    if(j>x2)&&(j<=x3) Xaxis(j)=(j-x2)*dx+X1+X2; end
    if(j>x3)&&(j<=x4) Xaxis(j)=(j-x3)*dx+X1+X2+X3; end
    if(j>x4)&&(j<=x5) Xaxis(j)=(j-x4)*dx+X1+X2+X3+X4; end
end
for i=1:NodeY
    if(i<=y1) Yaxis(i)=(i-1)*dy1; end
    if(i>y1)&&(i<=y2) Yaxis(i)=(i-y1)*dy2+Y1; end
    if(i>y2)&&(i<=y3) Yaxis(i)=(i-y2)*dy3+Y1+Y2; end
    if(i>y3)&&(i<=y4) Yaxis(i)=(i-y3)*dy4+Y1+Y2+Y3; end
end
% Plot model
PP(1,:)=[0, 0];
PP(2,:)=[X1, 0];
PP(3,:)=[X1+X2+X3+X4, 0];
PP(4,:)=[X1+X2+X3+X4+X5,0];
PP(5,:)=[X1, Y1];
PP(6,:)=[X1+X2+X3+X4, Y1];
PP(7,:)=[X1, Y1+Y2];
PP(8,:)=[X1+X2, Y1+Y2];
PP(9,:)=[X1+X2+X3, Y1+Y2];
PP(10,:)=[X1+X2+X3+X4, Y1+Y2];
PP(11,:)=[X1+X2, Y1+Y2+Y3];
PP(12,:)=[X1+X2+X3, Y1+Y2+Y3];
PP(13,:)=[0, Y1+Y2+Y3+Y4];
PP(14,:)=[X1, Y1+Y2+Y3+Y4];
PP(15,:)=[X1+X2+X3+X4, Y1+Y2+Y3+Y4];
PP(16,:)=[X1+X2+X3+X4+X5,Y1+Y2+Y3+Y4];
% Plot iron
IronX1=[PP(2,1) PP(5,1) PP(6,1) PP(3,1)];
IronY1=[PP(2,2) PP(5,2) PP(6,2) PP(3,2)];
IronX2=[PP(7,1),PP(8,1),PP(11,1),PP(12,1),PP(9,1),PP(10,1), PP(15,1),PP(14,1)];
IronY2=[PP(7,2),PP(8,2),PP(11,2),PP(12,2),PP(9,2),PP(10,2), PP(15,2),PP(14,2)];
% Plot copper
CopX1=[PP(8,1) PP(9,1) PP(12,1) PP(11,1)];
CopY1=[PP(8,2) PP(9,2) PP(12,2) PP(11,2)];
figure(1);
fill(IronX1,IronY1,'y'); hold on;
fill(IronX2,IronY2,'y'); hold on;

fill(CopX1,CopY1,'r');hold on;
axis equal; axis([0 7 0 6.5]);
contour(Xaxis,Yaxis,Az,30,'b');
disp('Amax='); disp(TotalFlux);
% Plot Meshes
figure(2);
fill(IronX1,IronY1,'y'); hold on;
fill(IronX2,IronY2,'y'); hold on;
fill(CopX1,CopY1,'r'); hold on;
axis equal; axis([0 7 0 6.5]);
for index=1:TotalMeshes
    i=MEM(index,1); j=MEM(index,2); k=MEM(index,3);
    xi=Nodes(i,1)/inTome; yi=Nodes(i,2)/inTome; xj=Nodes(j,1)/inTome;
    yj=Nodes(j,2)/inTome; xk=Nodes(k,1)/inTome; yk=Nodes(k,2)/inTome;
    line([xi,xj],[yi,yj]); line([xj,xk],[yj,yk]); line([xk,xi], [yk,yi]);
end
