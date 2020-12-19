%Ludovico Fossà 12/2020
%Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)

clear variables
close all
clc
options = optimset('Display','off');

%% INPUT
mach_exit=2.4;
gamma=1.4;
steps=50;
guess=0.375*pi/180;
plot_text=false;

%% OUTPUT VALUES
theta_wmax=prandtl_meyer(mach_exit,gamma)/2;
delta_theta=(theta_wmax-guess)/steps;
if(plot_text==true)
    col=sprintf('yellow');
else
    col=sprintf('black');
end

%number of points
nop=0;
for i=1:steps+1
    nop=nop+steps+2-(i-1);
end
data=zeros(nop,6); % point properties
XY=zeros(nop,2); %XY location - non dimensional throat radius=1
sel=@(x) steps+3+x*steps-(x+1).*x/2+2*x; %series

%% point properties and location
%FIRST ROUND: compute K- and K+ (characteristic constants)
for i=1:steps+1
    theta=guess+(i-1)*delta_theta;
    data(i,3)=theta;
    data(i,5)=fsolve(@(x) prandtl_meyer(x,gamma)-theta,1.1,options); %mach2
    data(i,6)=asin(1/data(i,5));
    data(i,4)=prandtl_meyer(data(i,5),gamma); %nu
    data(i,1)=data(i,3)+data(i,4); %K-
    data(i,2)=data(i,3)-data(i,4); %K+
end

%FIRST ROUND compute points on the first characteristic line
data(steps+2,:)=data(steps+1,:); %SIMPLE REGION
XY(1,1)=-1/tan(data(1,3)-data(1,6));
figure(1)
line([0 XY(1,1)],[1 XY(1,2)])
grid
hold on
for i=2:steps+1
    theta12=0.5*(data(i-1,3)+data(i,3));
    mu12=0.5*(data(i-1,6)+data(i,6));
    mu2=data(i,6);
    theta2=data(i,3);
    XY(i,1)=(XY(i-1,2)-1-tan(theta12+mu12)*XY(i-1,1))/(tan(theta2-mu2)-tan(theta12+mu12));
    XY(i,2)=XY(i-1,2)+tan(theta12+mu12)*(XY(i,1)-XY(i-1,1));
    line([XY(i-1,1) XY(i,1)],[XY(i-1,2) XY(i,2)])
    line([0 XY(i,1)],[1 XY(i,2)])
end
%first contour point
i=steps+2;
theta12=0.5*(data(i-1,3)+data(i,3));
mu12=0.5*(data(i-1,6)+data(i,6));
thetaai=0.5*(theta_wmax+data(i,3));
XY(i,1)=(XY(i-1,2)-1-tan(theta12+mu12)*XY(i-1,1))/(tan(thetaai)-tan(theta12+mu12));
XY(i,2)=XY(i-1,2)+tan(theta12+mu12)*(XY(i,1)-XY(i-1,1));
line([XY(i-1,1) XY(i,1)],[XY(i-1,2) XY(i,2)])
line([0 XY(i,1)],[1 XY(i,2)],'LineWidth',3,'Color',col) %contour

%SECOND ROUND: compute the non-simple region
for j=0:steps-1
   
   %CENTERLINE
   theta=0;
   data(sel(j),1)=data(j+2,1); %K-
   data(sel(j),4)=data(sel(j),1)+data(sel(j),3); %nu
   data(sel(j),2)=theta-data(sel(j),4); %K+
   data(sel(j),5)=fsolve(@(x) prandtl_meyer(x,gamma)-data(sel(j),4),1.1,options); %mach2
   data(sel(j),6)=asin(1./data(sel(j),5)); %mu
   
   %NON-SIMPLE REGION - moving along the K+ characteristic
   for i=1:steps-j-1
       data(sel(j)+i,3)=i*delta_theta; %no guess
       data(sel(j)+i,1)=data(j+i+2,1); %K-
       data(sel(j)+i,2)=data(sel(j),2); %K+
       data(sel(j)+i,4)=0.5*(data(sel(j)+i,1)-data(sel(j)+i,2)); %nu
       data(sel(j)+i,5)=fsolve(@(x) prandtl_meyer(x,gamma)-data(sel(j)+i,4),1.1,options); %mach2
       data(sel(j)+i,6)=asin(1./data(sel(j)+i,5)); %mu
   end
   data(sel(j)+steps-j,:)=data(sel(j)+steps-j-1,:); %SIMPLE REGION
   
    %DRAW CHARACTERISTIC LINES
    %draw centerline
    theta13=0.5*(data(sel(j-1)+1,3)+data(sel(j),3));
    mu13=0.5*(data(sel(j-1)+1,6)+data(sel(j),6));
    XY(sel(j),1)=XY(sel(j-1)+1,1)-XY(sel(j-1)+1,2)/tan(theta13-mu13);
    line([XY(sel(j-1)+1,1) XY(sel(j),1)],[XY(sel(j-1)+1,2) XY(sel(j),2)])
    for i=1:steps-j-1
        theta12=0.5*(data(sel(j)+i,3)+data(sel(j)+i-1,3));
        mu12=0.5*(data(sel(j)+i,6)+data(sel(j-1)+i-1,6));
        theta23=0.5*(data(sel(j)+i,3)+data(sel(j-1)+i+1,3));
        mu23=0.5*(data(sel(j)+i,6)+data(sel(j-1)+i+1,6));
        XY(sel(j)+i,1)=(XY(sel(j)+i-1,2)-XY(sel(j-1)+i+1,2)-tan(theta12+mu12)*XY(sel(j)+i-1,1)...
            +tan(theta23-mu23)*XY(sel(j-1)+i+1,1))/(tan(theta23-mu23)-tan(theta12+mu12));
        XY(sel(j)+i,2)=XY(sel(j)+i-1,2)+tan(theta12+mu12)*(XY(sel(j)+i,1)-XY(sel(j)+i-1,1));
        line([XY(sel(j-1)+1+i,1) XY(sel(j)+i,1)],[XY(sel(j-1)+i+1,2) XY(sel(j)+i,2)])
        line([XY(sel(j)+i-1,1) XY(sel(j)+i,1)],[XY(sel(j)+i-1,2) XY(sel(j)+i,2)])
    end
    %DRAW CHARACTERISTIC LINES - SIMPLE REGION - CONTOUR
    theta12=0.5*(data(sel(j+1)-1,3)+data(sel(j+1)-2,3));
    mu12=0.5*(data(sel(j+1)-1,6)+data(sel(j+1)-2,6));
    theta23=0.5*(data(sel(j)-1,3)+data(sel(j+1)-1,3));
    XY(sel(j+1)-1,1)=(XY(sel(j+1)-2,2)-XY(sel(j)-1,2)-tan(theta12+mu12)*XY(sel(j+1)-2,1)...
        +tan(theta23)*XY(sel(j)-1,1))/(tan(theta23)-tan(theta12+mu12));
    XY(sel(j+1)-1,2)=XY(sel(j+1)-2,2)+tan(theta12+mu12)*(XY(sel(j+1)-1,1)-XY(sel(j+1)-2,1));
    line([XY(sel(j+1)-2,1) XY(sel(j+1)-1,1)],[XY(sel(j+1)-2,2) XY(sel(j+1)-1,2)])
    line([XY(sel(j)-1,1) XY(sel(j+1)-1,1)],[XY(sel(j)-1,2) XY(sel(j+1)-1,2)],'LineWidth',3,'Color',col)%CONTOUR
end
daspect([1 1 0.2778])
title('Minimum-length nozzle contour','FontSize',20,'FontWeight','bold')
xlabel('x/r_{throat}','FontSize',16,'FontWeight','bold')
ylabel('y/r_{throat}','FontSize',16,'FontWeight','bold')

%% EXPORT
export=zeros(steps+2,3);
export(:,1)=1;
export(1,:)=[1 0 1];
for j=-1:steps-1
    export(j+3,2)=XY(sel(j+1)-1,1);
    export(j+3,3)=XY(sel(j+1)-1,2);
end

format short
writematrix(export,'nozzle.txt','Delimiter','tab')

%% ADD TEXT TO THE PLOT
figure(1)
if(plot_text==true)
    for i=1:size(XY,1)
        stringa=sprintf('%d',i);
        text(XY(i,1),XY(i,2),stringa)
    end
end
