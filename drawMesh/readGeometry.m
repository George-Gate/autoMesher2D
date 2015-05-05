fileID=fopen('../geometry.txt');
ng=fscanf(fileID,'n=%d\n',1);
fgetl(fileID);
g=fscanf(fileID,'%g\n',[3,ng+1]);
g=g';

plot(g(:,1),g(:,2),g(g(:,3)>0,1),g(g(:,3)>0,2),'o');