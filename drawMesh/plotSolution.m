%% read data
fileID=fopen('../solution.txt');
n=fscanf(fileID,'N=%d\n',1);
fgetl(fileID); % skip title
xyu=fscanf(fileID,'%g\n',[3,n]);
xyu=xyu';

%% interpolation
u=scatteredInterpolant(xyu(:,1),xyu(:,2),xyu(:,3),'linear','none');
xmin=min(xyu(:,1));xmax=max(xyu(:,1));
ymin=min(xyu(:,2));ymax=max(xyu(:,2));
[x,y]=ndgrid(linspace(xmin,xmax),linspace(ymin,ymax));

%% plot
readGeometry;
hold on;
contour(x,y,u(x,y),10);
hold off;