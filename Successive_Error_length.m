clc
clear
close all
%input Datas
Re=input('Reynolds number= ');
t=input('time= ');

W=20;
z=zeros(1,W-1);
Z=zeros(1,W-1);

for w=0:W
    if w>1
       ...Midpoint temperature: 
       M=V(2,0.5*(J+2));
       z(1,w-1)=h;
    end
...Number of pieces created on the length [0,1]
J=10+(2*w+2);
h=1/J;
...Time step
ta=10^-3/2;
t=0:ta:t;
... Re=U*L/nu
U=1;
nu=U*1/Re;
...Velocity field u_n
    for m=1:max(size(t))
    if m==1
u_n=zeros(J+1);
u_n=[ones(1,J+1); u_n(2:end,1:end)];
...Velocity field v_n
v_n=zeros(J+1);
...Velocity field u*
    else
    u_n=u_new;
    v_n=v_new;
    end
u_star=zeros(J+1);
v_star=zeros(J+1);
for i=2:J
  for j=2:J
      u_star(i,j)=-(ta/(2*h))*(u_n(i,j+1)^2-u_n(i,j)^2+2*u_n(i-1,j)*v_n(i-1,j)-2*u_n(i,j)*v_n(i,j))+(nu*ta/h^2)*(u_n(i,j-1)-4*u_n(i,j)+u_n(i,j+1)+u_n(i+1,j)+u_n(i-1,j))+u_n(i,j);
      v_star(i,j)=-(ta/(2*h))*(2*u_n(i,j+1)*v_n(i,j+1)-2*u_n(i,j)*v_n(i,j)+v_n(i-1,j)^2-v_n(i,j)^2)+(nu*ta/h^2)*(v_n(i,j-1)-4*v_n(i,j)+v_n(i,j+1)+v_n(i+1,j)+v_n(i-1,j))+v_n(i,j);
  end
end
u_star=(u_star(2:end-1,2:end-1));
v_star=(v_star(2:end-1,2:end-1));
...Create a coefficient matrix A
A=zeros((J-1)^2);
for i=1:(J-1)^2
    for j=1:(J-1)^2
        if i==j
            A(i,j)=4;
        elseif  i==j+1
            if mod(j,J-1)==0
               A(i,j)=-0;
            else
               A(i,j)=-1;
            end
        elseif i==j-1
            if mod(i,J-1)==0
               A(i,j)=0;
            else
               A(i,j)=-1;
            end
        elseif i-j==J-1
            A(i,j)=-1;
        elseif j-i==J-1
            A(i,j)=-1;        
        end
    end
end
...Create a coefficient matrix b
g=zeros(J-1);
for j=1:J-2
  for i=2:J-1
   g(i,j)=(1/(ta*h))*(u_star(i,j+1)-u_star(i,j)+v_star(i-1,j)-v_star(i,j));
  end
end
g(1,1:J-2)=2.*g(2,1:J-2)-g(3,1:J-2);
g(1:J-1,J-1)=2.*g(1:J-1,J-2)-g(1:J-1,J-3);
g=reshape(rot90(g,3),(J-1)^2,1);
b=-h^2.*g;
...P' calculation
P=linsolve(A,b);
P=rot90(reshape(P,J-1,J-1));
...u_n+1 & v_n+1 calculation
[px, py]=gradient(P);
u_new=u_star-(ta.*px);
u_new=[ones(1,J+1) ;zeros(J-1,1) u_new zeros(J-1,1); zeros(1,J+1)];
v_new=v_star-(ta.*py);
v_new=[zeros(1,J+1) ;zeros(J-1,1) v_new zeros(J-1,1); zeros(1,J+1)];
V=sqrt(u_new.^2+v_new.^2);
    end
...Calculation of successive error
if w>1
Z(1,w-1)=abs(V(2,0.5*(J+2))-M);
end
end
figure(1)
loglog(z,Z,'o')
xlabel('h')
ylabel('Error')
figure(2)
slope=gradient(log(Z))./gradient(log(z));
plot(z,slope,'o')
xlabel('h')
ylabel('Slope')