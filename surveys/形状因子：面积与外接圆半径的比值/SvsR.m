%%
angle=0.01:0.05:pi-0.01;
a=1;
b=0.01:0.01:10;

for the=angle
    S=0.5*b*sin(the);
    c=sqrt(a^2+b.^2-2*a*b*cos(the));
    R=a*b.*c./(4*S);
    plot(b,S./(pi*R.^2));
    title(['\theta=',num2str(round(the/pi*180)),'°']);
    set(gca,'Ylim',[0 0.5]);
    pause(0.2);
end

%% (init) S./R^2 的角度分布
[the1,the2]=meshgrid(0.01:0.005:pi-0.01,0.01:0.005:pi-0.01);
a=1;
h=a./(cot(the1)+cot(the2));
S=0.5*a*h;
b=h./sin(the1);
c=h./sin(the2);
R=a.*b.*c./(4*S);
R(the1+the2>=pi)=NaN;
mesh(the1/pi*180,the2/pi*180,S./(pi*R.^2));
xlabel('\theta_1');
ylabel('\theta_2');
zlabel('S/\piR^2');

%% S./R^2 vs max angle
maxAngle=max(max(the1,the2),pi-the1-the2);
ratio=S./(pi*R.^2);
angleList=pi/3:0.005:pi;
ub=angleList;
lb=angleList;
for i=1:length(angleList)
    sample=ratio(abs(maxAngle-angleList(i))<0.0025);
    if (isempty(sample))
        ub(i)=NaN;
        lb(i)=NaN;
    else
        ub(i)=max(sample);
        lb(i)=min(sample);
    end
end

plot(angleList/pi*180,ub,angleList/pi*180,lb);
xlabel('Max Angle');
ylabel('S/\piR^2 Area');

%% S./R^2 vs min angle
minAngle=min(min(the1,the2),pi-the1-the2);
ratio=S./(pi*R.^2);
angleList=0:0.005:pi/3;
ub=angleList;
lb=angleList;
for i=1:length(angleList)
    sample=ratio(abs(minAngle-angleList(i))<0.0025);
    if (isempty(sample))
        ub(i)=NaN;
        lb(i)=NaN;
    else
        ub(i)=max(sample);
        lb(i)=min(sample);
    end
end

plot(angleList/pi*180,ub,angleList/pi*180,lb);
xlabel('Min Angle');
ylabel('S/\piR^2 Area');

%% maxEdge/R   minEdge/R
minEdge=min(min(a,b),c);
maxEdge=max(max(a,b),c);
mesh(the1/pi*180,the2/pi*180,maxEdge./R);
xlabel('\theta_1');
ylabel('\theta_2');
zlabel('R./minEdge');
