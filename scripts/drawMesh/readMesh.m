readGeometry;
%% read data
fileID=fopen('../mesh.txt');
mesh=struct([]);
i=1;
while (1)
    tag=fgetl(fileID);
    if (feof(fileID)) 
        break;
    end
    mesh(i).tag=deblank(tag);
    % read vertex
    nv=fscanf(fileID,'N_Vertex=%d\n',1);
    fgetl(fileID); % skip title
    vertex=fscanf(fileID,'%e %e %d',[3,nv]);
    mesh(i).vertex=vertex';
    mesh(i).nv=nv;
       
    % read edge
    fgetl(fileID);fgetl(fileID);   % skip blank line
    ne=fscanf(fileID,'N_Edge=%d\n',1);
    fgetl(fileID);  % skip title
    edge=fscanf(fileID,'%d %d',[2,ne]);
    mesh(i).edge=edge';
    mesh(i).ne=ne;
    
    % read surface
    fgetl(fileID);fgetl(fileID);   % skip blank line
    ns=fscanf(fileID,'N_Surface=%d\n',1);
    fgetl(fileID);  % skip title
    surface=fscanf(fileID,'%d %d %d %d %d %d',[6,ns]);
    mesh(i).surface=surface';
    mesh(i).ns=ns;
    
    % read boundary
    fgetl(fileID);fgetl(fileID);   % skip blank line
    nb=fscanf(fileID,'N_Boundary=%d\n',1);
    fgetl(fileID);  % skip title
    boundary=fscanf(fileID,'%d',nb);
    mesh(i).boundary=boundary;
    mesh(i).nb=nb;
    
    fgetl(fileID);fgetl(fileID);  % skip blank line
    i=i+1;
end
fclose(fileID);

%% plot
for i=1:length(mesh)
   ne=mesh(i).ne;
   hold off;
   plot(g(:,1),g(:,2));
   hold on;
   scatter(g(g(:,3)>0,1),g(g(:,3)>0,2),'red');
   for j=1:ne
       line(mesh(i).vertex(mesh(i).edge(j,1:2),1),mesh(i).vertex(mesh(i).edge(j,1:2),2));
   end
   scatter(mesh(i).vertex(:,1),mesh(i).vertex(:,2),'.','blue');
   title(mesh(i).tag);
   pause(0.1);
end
