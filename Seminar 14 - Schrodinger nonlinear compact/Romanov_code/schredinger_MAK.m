function [U1,V1, ep1, ep2, integral, integral1,integral2,  test1] = schredinger_MAK(h, D, t, y, x)

% y -- true solution (from rynge.m)
% x -- true solution (from rynge.m)

% n1 = 200;
% m1 = 1024;
n1 = ceil(2*pi/h);
m1 = ceil(2*pi/t);

%X = linspace(-10, 10, n1+1);

%X = linspace(0, 2*pi, n1+1);
X = x;
%X = linspace(0, 1, n1+1);
X = X(1:length(x)-1)';


u0 = cos(2*pi/(2.3)*X).*y(1:length(y)-1, 1)'+0.3;
v0 = sin(2*pi/(2.3)*X).*y(1:length(y)-1, 1)'+0.3;

u0 = zeros(length(x)-1, 1)+1;
v0 = zeros(length(x)-1, 1);

u01 = zeros(length(x)-1, 1);
v01 = zeros(length(x)-1, 1);

n1 = length(x)-1;
integral = zeros(4*m1+1, 1);

integral1 = zeros(4*m1+1, 1);
integral2 = zeros(4*m1+1, 1);

u01(1) = (u0(2) - u0(n1))/2;
v01(1) = (v0(2) - v0(n1))/2;

for i = 2:length(u0)-1
    u01(i) = (u0(i+1)-u0(i-1))/2;
    v01(i) = (v0(i+1)-v0(i-1))/2;
end

u01(n1) = (u0(1) - u0(n1-1))/2;
v01(n1) = (v0(1) - v0(n1-1))/2;

cur1 = -u0.*v01+v0.*u01;

cur = u0.^2+v0.^2;
cur2 = (u0.^2+v0.^2).^2+u01.^2+v01.^2;

summ = 0;
summ1 = 0;
summ2 = 0;
for i = 1:length(cur)
    summ = summ + cur(i);
    summ1 = summ1 + cur1(i);
    summ2 = summ2 + cur2(i);
end
integral(1) = summ;
integral1(1) = summ1;
integral2(1) = summ2;

F = zeros(n1, 4*m1+1);


%ep1 = zeros(n1, m1+1);
ep1 = zeros(n1, 4*m1+1);
ep2 = zeros(n1, 4*m1+1);
%U_test(:, 1) = u0;
U1 = zeros(n1, 4*m1+1);
U1(:, 1) = u0;

V1 = zeros(n1, 4*m1+1);
V1(:, 1) = v0;


 for i = 1:4*m1
     u = zeros(n1, 1);
     v = zeros(n1, 1);
     u_dop = zeros(n1, 1);
     f = zeros(n1, 1);
     v_dop = zeros(n1, 1);

     for j = 1:n1
     
           if j == 1
                 
                  
                  u_dop0= u0(n1)-D*t/(h^2)*(v0(j)-2*v0(n1)+v0(n1-1))-(u0(n1)^2+v0(n1)^2)*v0(n1)*t;
                  u_dop(j)=u0(j)-D*t/(h^2)*(v0(j+1)-2*v0(j)+v0(n1))-(u0(j)^2+v0(j)^2)*v0(j)*t;
                  u_dop2 = u0(j+1)-D*t/(h^2)*(v0(j+2)-2*v0(j+1)+v0(j))-(u0(j+1)^2+v0(j+1)^2)*v0(j+1)*t;
                  
                  v_dop0= v0(n1)+D*t/(h^2)*(u0(j)-2*u0(n1)+u0(n1-1))+(u0(n1)^2+v0(n1)^2)*u0(n1)*t;
                  v_dop(j)=v0(j)+D*t/(h^2)*(u0(j+1)-2*u0(j)+u0(n1))+(u0(j)^2+v0(j)^2)*u0(j)*t;
                  v_dop2 = v0(j+1)+D*t/(h^2)*(u0(j+2)-2*u0(j+1)+u0(j))+(u0(j+1)^2+v0(j+1)^2)*u0(j+1)*t;
                  
                  
                  u(j) = u0(j)/2+u_dop(j)/2-D*t/(2*h^2)*(v_dop2-2*v_dop(j)+v_dop0)-(u_dop(j)^2+v_dop(j)^2)*v_dop(j)*t/2;
                  v(j) = v0(j)/2+v_dop(j)/2+D*t/(2*h^2)*(u_dop2-2*u_dop(j)+u_dop0)+(u_dop(j)^2+v_dop(j)^2)*u_dop(j)*t/2;

           
           elseif j == n1
               
               
                  u_dop0= u0(j-1)-D*t/(h^2)*(v0(j)-2*v0(j-1)+v0(j-2))-(u0(j-1)^2+v0(j-1)^2)*v0(j-1)*t;
                  u_dop(j)=u0(j)-D*t/(h^2)*(v0(1)-2*v0(j)+v0(j-1))-(u0(j)^2+v0(j)^2)*v0(j)*t;
                  u_dop2 = u0(1)-D*t/(h^2)*(v0(2)-2*v0(1)+v0(j))-(u0(1)^2+v0(1)^2)*v0(1)*t;
                  
                  v_dop0= v0(j-1)+D*t/(h^2)*(u0(j)-2*u0(j-1)+u0(j-2))+(u0(j-1)^2+v0(j-1)^2)*u0(j-1)*t;
                  v_dop(j)=v0(j)+D*t/(h^2)*(u0(1)-2*u0(j)+u0(j-1))+(u0(j)^2+v0(j)^2)*u0(j)*t;
                  v_dop2 = v0(1)+D*t/(h^2)*(u0(2)-2*u0(1)+u0(j))+(u0(1)^2+v0(1)^2)*u0(1)*t;

                  
                  
                  u(j) = u0(j)/2+u_dop(j)/2-D*t/(2*h^2)*(v_dop2-2*v_dop(j)+v_dop0)-(u_dop(j)^2+v_dop(j)^2)*v_dop(j)*t/2;
                  v(j) = v0(j)/2+v_dop(j)/2+D*t/(2*h^2)*(u_dop2-2*u_dop(j)+u_dop0)+(u_dop(j)^2+v_dop(j)^2)*u_dop(j)*t/2;
               
           else

               u_dop(j)=u0(j)-D*t/(h^2)*(v0(j+1)-2*v0(j)+v0(j-1))-(u0(j)^2+v0(j)^2)*v0(j)*t;
               v_dop(j)=v0(j)+D*t/(h^2)*(u0(j+1)-2*u0(j)+u0(j-1))+(u0(j)^2+v0(j)^2)*u0(j)*t;
               if j == 2

                  u_dop0= u0(j-1)-D*t/(h^2)*(v0(j)-2*v0(j-1)+v0(n1))-(u0(j-1)^2+v0(j-1)^2)*v0(j-1)*t;
                  u_dop2 = u0(j+1)-D*t/(h^2)*(v0(j+2)-2*v0(j+1)+v0(j))-(u0(j+1)^2+v0(j+1)^2)*v0(j+1)*t;
                  
                  v_dop0= v0(j-1)+D*t/(h^2)*(u0(j)-2*u0(j-1)+u0(n1))+(u0(j-1)^2+v0(j-1)^2)*u0(j-1)*t;
                  v_dop2 = v0(j+1)+D*t/(h^2)*(u0(j+2)-2*u0(j+1)+u0(j))+(u0(j+1)^2+v0(j+1)^2)*u0(j+1)*t;

               elseif j== n1-1

                  u_dop0= u0(j-1)-D*t/(h^2)*(v0(j)-2*v0(j-1)+v0(j-2))-(u0(j-1)^2+v0(j-1)^2)*v0(j-1)*t;
                  u_dop2 = u0(j+1)-D*t/(h^2)*(v0(1)-2*v0(j+1)+v0(j))-(u0(j+1)^2+v0(j+1)^2)*v0(j+1)*t;
                  
                  v_dop0= v0(j-1)+D*t/(h^2)*(u0(j)-2*u0(j-1)+u0(j-2))+(u0(j-1)^2+v0(j-1)^2)*u0(j-1)*t;
                  v_dop2 = v0(j+1)+D*t/(h^2)*(u0(1)-2*u0(j+1)+u0(j))+(u0(j+1)^2+v0(j+1)^2)*u0(j+1)*t;

               else
              

                  u_dop0= u0(j-1)-D*t/(h^2)*(v0(j)-2*v0(j-1)+v0(j-2))-(u0(j-1)^2+v0(j-1)^2)*v0(j-1)*t;
                  u_dop2 = u0(j+1)-D*t/(h^2)*(v0(j+2)-2*v0(j+1)+v0(j))-(u0(j+1)^2+v0(j+1)^2)*v0(j+1)*t;
                  
                  v_dop0= v0(j-1)+D*t/(h^2)*(u0(j)-2*u0(j-1)+u0(j-2))+(u0(j-1)^2+v0(j-1)^2)*u0(j-1)*t;
                  v_dop2 = v0(j+1)+D*t/(h^2)*(u0(j+2)-2*u0(j+1)+u0(j))+(u0(j+1)^2+v0(j+1)^2)*u0(j+1)*t;
                    

               end
               
                  u(j) = u0(j)/2+u_dop(j)/2-D*t/(2*h^2)*(v_dop2-2*v_dop(j)+v_dop0)-(u_dop(j)^2+v_dop(j)^2)*v_dop(j)*t/2;
                  v(j) = v0(j)/2+v_dop(j)/2+D*t/(2*h^2)*(u_dop2-2*u_dop(j)+u_dop0)+(u_dop(j)^2+v_dop(j)^2)*u_dop(j)*t/2;
               
           end   
        
            
     end

     
      u0 = u;
      v0 = v;
      
      u01 = zeros(length(x)-1, 1);
      v01 = zeros(length(x)-1, 1);
      
      u01(1) = (u0(2) - u0(n1))/2;
      v01(1) = (v0(2) - v0(n1))/2;
      
      for j = 2:length(u0)-1
        u01(j) = (u0(j+1)-u0(j-1))/2;
        v01(j) = (v0(j+1)-v0(j-1))/2;
      end
      
      u01(n1) = (u0(1) - u0(n1-1))/2;
      v01(n1) = (v0(1) - v0(n1-1))/2;
      
      E1 = diag(ones(n1, 1)*2/3)+diag(ones(n1-1,1)*1/6,1)+diag(ones(n1-1,1)*1/6,-1);
      E1(1, n1) = 1/6;
      E1(n1, 1) = 1/6;
      ans_cur_u = inv(E1)*u01;
      ans_cur_v = inv(E1)*v01;
      cur = u0.^2+v0.^2;
           
      cur1 = -u0.*ans_cur_v+v0.*ans_cur_u;
      cur2 = (u0.^2+v0.^2).^2+u01.^2+v01.^2;
      summ = 0;
      summ1 = 0;
      summ2 = 0;
      
      for j = 1:length(cur)
            summ = summ + cur(j);
            summ1 = summ1 + cur1(j);
            summ2 = summ2 + cur2(j);
      end

      integral(i+1) = summ;
      integral1(i+1) = summ1;
      integral2(i+1) = summ2;

      F(:, i+1) = f;
      if i == 2
          test1 = E1;
      end
          
      U1(:, i+1) = u0;
      V1(:, i+1) = v0;
      

   
    
 end
