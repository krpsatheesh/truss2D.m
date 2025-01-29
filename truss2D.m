% getting nodes
noofnodes= input('Enter number of nodes: ');
%nodedata: nodenum,x,y,xrestraint,yrestraint,fx,fy
nodedata(noofnodes,7)=[0];
for temp=1:noofnodes
 nodedata(temp,1)=temp;
 fprintf('Enter the details of node %d, like:',temp);
 nodedata(temp,2:7)=input('[x,y,xrestraint,yrestraint,fx,fy]: ');
end

%getting members
noofmembers= input('Enter number of members: ');
%memberdata: memebernum,i,j,A,E
memberdata(noofmembers,5)=[0];

for temp=1:noofmembers
 memberdata(temp,1)=temp;
 fprintf('Enter the details of member %d, like:',temp);
 memberdata(temp,2:5)=input('[inode,jnode,A,E]: ');
end

%calculating length of member
%memberdata: memebernum,i,j,A,E,L
memberdata(:,6)=[0];
for temp=1:noofmembers
  xi=nodedata(memberdata(temp,2),2);
  xj=nodedata(memberdata(temp,3),2);
  yi=nodedata(memberdata(temp,2),3);
  yj=nodedata(memberdata(temp,3),3);
  memberdata(temp,6)=sqrt((xj-xi)^2+(yj-yi)^2);
end

%calculating direction cosines
%memberdata: memebernum,i,j,A,E,L,c,s
memberdata(:,7:8)=[0];
for temp=1:noofmembers
  xi=nodedata(memberdata(temp,2),2);
  xj=nodedata(memberdata(temp,3),2);
  yi=nodedata(memberdata(temp,2),3);
  yj=nodedata(memberdata(temp,3),3);
  memberdata(temp,7)=(xj-xi)/memberdata(temp,6);
  memberdata(temp,8)=(yj-yi)/memberdata(temp,6);
end

%assembling global stiffness matrix
globalk(noofnodes*2,noofnodes*2)=[0];
for temp=1:noofmembers
  i=memberdata(temp,2);
  j=memberdata(temp,3);
  A=memberdata(temp,4);
  E=memberdata(temp,5);
  L=memberdata(temp,6);
  c2=memberdata(temp,7)^2;
  s2=memberdata(temp,8)^2;
  cs=memberdata(temp,7)*memberdata(temp,8);
  globalk(2*i-1:2*i,2*i-1:2*i)=globalk(2*i-1:2*i,2*i-1:2*i)+A*(E/L)*[
  c2 cs
  cs s2];
  globalk(2*i-1:2*i,2*j-1:2*j)=globalk(2*i-1:2*i,2*j-1:2*j)+A*(E/L)*[
  -c2 -cs
  -cs -s2];
  globalk(2*j-1:2*j,2*i-1:2*i)=globalk(2*j-1:2*j,2*i-1:2*i)+A*(E/L)*[
  -c2 -cs
  -cs -s2];
  globalk(2*j-1:2*j,2*j-1:2*j)=globalk(2*j-1:2*j,2*j-1:2*j)+A*(E/L)*[
  c2 cs
  cs s2];
end

%initialize u and f vectors
u(noofnodes*2)=[0];
f(noofnodes*2)=[0];
u=transpose(u);
f=transpose(f);
for temp=1:noofnodes
 f(temp*2-1)=nodedata(temp,6);
 f(temp*2)=nodedata(temp,7);
end

%matrix reduction
reducedk=globalk;  %initializing reducedk

temp2=0;
for temp=1:noofnodes
  if(nodedata(temp,4)==1)  %Identifying restrained Coordinates
    temp2=temp2+1;
    zeroref(temp2)=temp*2-1;
  end
   if(nodedata(temp,5)==1)
    temp2=temp2+1;
    zeroref(temp2)=temp*2;
   end
end

reducedk(zeroref,:)=[];
reducedk(:,zeroref)=[];

fprintf('Global Stiffness Matrix is:\n');
disp(globalk);
fprintf('Reduced Stiffness Matrix is:\n');
disp(reducedk);

%solving for unknown displacements
knownf=f;
knownf(zeroref,:)=[];
unknownu=inv(reducedk)*knownf;

%complete displacement vector
temp2=1;
for temp=1:noofnodes*2
  if (any(zeroref(:)==temp))
    u(temp)=0;
  else
    u(temp)=unknownu(temp2);
    temp2=temp2+1;
  end
end

%display deflection vector
disp("unknown displacements are:");
disp(unknownu);

%solving for unknown reactions
globalf=globalk*u;

%display reaction vector
disp("unknown reactions are:");
disp(globalf(zeroref));

%member forces calculation
for temp=1:noofmembers
  i=memberdata(temp,2);
  j=memberdata(temp,3);
  A=memberdata(temp,4);
  E=memberdata(temp,5);
  L=memberdata(temp,6);
  c=memberdata(temp,7);
  s=memberdata(temp,8);
  ui=u(i*2-1);
  vi=u(i*2);
  uj=u(j*2-1);
  vj=u(j*2);
  memberforces(temp)= A*(E/L)*[c s -c -s]*transpose([ui vi uj vj]);
end

memberforces=transpose(memberforces);
disp("Force in members are: ");
disp(memberforces);
