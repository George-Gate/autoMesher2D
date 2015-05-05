%% read data
readGeometry;
fileID=fopen('../finalMesh.txt');
finalMesh=struct();
% read vertex
nv=fscanf(fileID,'N_Vertex=%d\n',1);
fgetl(fileID); % skip title
vertex=fscanf(fileID,'%e %e %d',[3,nv]);
finalMesh.vertex=vertex';
finalMesh.nv=nv;

% read edge
fgetl(fileID);fgetl(fileID);   % skip blank line
ne=fscanf(fileID,'N_Edge=%d\n',1);
fgetl(fileID);  % skip title
edge=fscanf(fileID,'%d %d',[2,ne]);
finalMesh.edge=edge';
finalMesh.ne=ne;

% read surface
fgetl(fileID);fgetl(fileID);   % skip blank line
ns=fscanf(fileID,'N_Surface=%d\n',1);
fgetl(fileID);  % skip title
surface=fscanf(fileID,'%d %d %d %d %d %d',[6,ns]);
finalMesh.surface=surface';
finalMesh.ns=ns;

% read boundary
fgetl(fileID);fgetl(fileID);   % skip blank line
nb=fscanf(fileID,'N_Boundary=%d\n',1);
fgetl(fileID);  % skip title
boundary=fscanf(fileID,'%d',nb);
finalMesh.boundary=boundary;
finalMesh.nb=nb;

fgetl(fileID);fgetl(fileID);  % skip blank line

fclose(fileID);

%% plot
ne=finalMesh.ne;
hold off;
plot(g(:,1),g(:,2));
hold on;
scatter(g(g(:,3)>0,1),g(g(:,3)>0,2),'red');
for j=1:ne
   line(finalMesh.vertex(finalMesh.edge(j,1:2),1),finalMesh.vertex(finalMesh.edge(j,1:2),2));
end
scatter(finalMesh.vertex(:,1),finalMesh.vertex(:,2),'.','blue');
% plot vertex id
nv=finalMesh.nv;
for j=1:nv
    text(finalMesh.vertex(j,1),finalMesh.vertex(j,2),num2str(finalMesh.vertex(j,3)),'VerticalAlignment','Bottom');
end
axis off
