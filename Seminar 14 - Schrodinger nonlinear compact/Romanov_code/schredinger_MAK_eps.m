function [U1,V1, ep1, ep2, integral, integral1, integral2, anss, spectr] = schredinger_MAK_eps(h, D, t, y, x)

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

integral = zeros(2*m1+1, 1);
integral1 = zeros(4*m1+1, 1);
integral2 = zeros(4*m1+1, 1);
spectr = zeros(4*m1+1, 1);
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
     f11= zeros(2*n1, 1);
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
     
     a = -0.51;
     
     a1 = a;
     b1 = 1;
     c1 = a;
     a0 = -a;
     b0 = -1;
     c0 = -a;
     f1 = (D*t+2*D*a*t)/(2*h^2);
     e1 = -(D*t+2*D*a*t)/(h^2);
     d1 = (D*t+2*D*a*t)/(2*h^2);
     
     f0 = (D*t+2*D*a*t)/(2*h^2);
     e0 = -(D*t+2*D*a*t)/(h^2);
     d0 = (D*t+2*D*a*t)/(2*h^2);
     p1 = -a*t/2;
     q1 = -t/2;
     r1 = -a*t/2;
     p0 = -a*t/2;
     q0 = -t/2;
     r0 = -a*t/2;
     
     c = -0.51;
%      A1 = -(D*t + 2*D*c*t)/(2*h^2);
%      B1 = (D*t + 2*D*c*t)/h^2;
%      C1 = -(D*t + 2*D*c*t)/(2*h^2);
%      A0 = -(D*t + 2*D*c*t)/(2*h^2);
%      B0 = (D*t + 2*D*c*t)/h^2;
%      C0 = -(D*t + 2*D*c*t)/(2*h^2);

     A1 = c;
     B1 = 1;
     C1 = c;
     A0 = -c;
     B0 = -1;
     C0 = -c;
     
%      F1 = c;
%      E1 = 1;
%      D1 =c;
%      F0 = -c;
%      E0 = -1;
%      D0 = -c;

     F1 = -(D*t + 2*D*c*t)/(2*h^2);
     E1 = (D*t + 2*D*c*t)/h^2;
     D1 =-(D*t + 2*D*c*t)/(2*h^2);
     F0 = -(D*t + 2*D*c*t)/(2*h^2);
     E0 = (D*t + 2*D*c*t)/h^2;
     D0 = -(D*t + 2*D*c*t)/(2*h^2);
     
     P1 = c*t/2;
     Q1 = t/2;
     R1 = c*t/2;
     P0 = c*t/2;
     Q0 = t/2;
     R0 = c*t/2;


      

      for j = 1:n1
        
            if j ==1 
                
                f11(j) = -a1*u(n1)-b1*u(j)-c1*u(j+1)-a0*u0(n1)-b0*u0(j)-c0*u0(j+1)-d1*v(n1)-e1*v(j)-f1*v(j+1)-d0*v0(n1)-e0*v0(j)-f0*v0(j+1)+p1*(u(n1)^2*v(n1)+v(n1)^3)+q1*(u(j)^2*v(j)+v(j)^3)+r1*(u(j+1)^2*v(j+1)+v(j+1)^3)+p0*(u0(n1)^2+v0(n1)^2)*v0(n1)+q0*(u0(j)^2+v0(j)^2)*v0(j)+r0*(u0(j+1)^2+v0(j+1)^2)*v0(j+1);
                f11(n1+j) = -A1*v(n1)-B1*v(j)-C1*v(j+1)-A0*v0(n1)-B0*v0(j)-C0*v0(j+1)-D1*u(n1)-E1*u(j)-F1*u(j+1)-D0*u0(n1)-E0*u0(j)-F0*u0(j+1)+P1*(v(n1)^2*u(n1)+u(n1)^3)+Q1*(v(j)^2*u(j)+u(j)^3)+R1*(v(j+1)^2*u(j+1)+u(j+1)^3)+P0*(u0(n1)^2+v0(n1)^2)*u0(n1)+Q0*(u0(j)^2+v0(j)^2)*u0(j)+R0*(u0(j+1)^2+v0(j+1)^2)*u0(j+1);

            elseif j == n1
           
                f11(j) = -a1*u(j-1)-b1*u(j)-c1*u(1)-a0*u0(j-1)-b0*u0(j)-c0*u0(1)-d1*v(j-1)-e1*v(j)-f1*v(1)-d0*v0(j-1)-e0*v0(j)-f0*v0(1)+p1*(u(j-1)^2*v(j-1)+v(j-1)^3)+q1*(u(j)^2*v(j)+v(j)^3)+r1*(u(1)^2*v(1)+v(1)^3)+p0*(u0(j-1)^2+v0(j-1)^2)*v0(j-1)+q0*(u0(j)^2+v0(j)^2)*v0(j)+r0*(u0(1)^2+v0(1)^2)*v0(1);
                f11(n1+j) = -A1*v(j-1)-B1*v(j)-C1*v(1)-A0*v0(j-1)-B0*v0(j)-C0*v0(1)-D1*u(j-1)-E1*u(j)-F1*u(1)-D0*u0(j-1)-E0*u0(j)-F0*u0(1)+P1*(v(j-1)^2*u(j-1)+u(j-1)^3)+Q1*(v(j)^2*u(j)+u(j)^3)+R1*(v(1)^2*u(1)+u(1)^3)+P0*(u0(j-1)^2+v0(j-1)^2)*u0(j-1)+Q0*(u0(j)^2+v0(j)^2)*u0(j)+R0*(u0(1)^2+v0(1)^2)*u0(1);
                
            else

                f11(j) = -a1*u(j-1)-b1*u(j)-c1*u(j+1)-a0*u0(j-1)-b0*u0(j)-c0*u0(j+1)-d1*v(j-1)-e1*v(j)-f1*v(j+1)-d0*v0(j-1)-e0*v0(j)-f0*v0(j+1)+p1*(u(j-1)^2*v(j-1)+v(j-1)^3)+q1*(u(j)^2*v(j)+v(j)^3)+r1*(u(j+1)^2*v(j+1)+v(j+1)^3)+p0*(u0(j-1)^2+v0(j-1)^2)*v0(j-1)+q0*(u0(j)^2+v0(j)^2)*v0(j)+r0*(u0(j+1)^2+v0(j+1)^2)*v0(j+1);
                f11(n1+j) = -A1*v(j-1)-B1*v(j)-C1*v(j+1)-A0*v0(j-1)-B0*v0(j)-C0*v0(j+1)-D1*u(j-1)-E1*u(j)-F1*u(j+1)-D0*u0(j-1)-E0*u0(j)-F0*u0(j+1)+P1*(v(j-1)^2*u(j-1)+u(j-1)^3)+Q1*(v(j)^2*u(j)+u(j)^3)+R1*(v(j+1)^2*u(j+1)+u(j+1)^3)+P0*(u0(j-1)^2+v0(j-1)^2)*u0(j-1)+Q0*(u0(j)^2+v0(j)^2)*u0(j)+R0*(u0(j+1)^2+v0(j+1)^2)*u0(j+1);

            end
      end
    
    E = diag(ones(2*n1, 1)*0);
    
    %first square
    for j = 1:n1-1
        
        E(j+1, j) = a1-2*p1*u(j)*v(j);

     end

     for j = 2:n1
         
        E(j-1, j) = c1-2*r1*u(j)*v(j);

     end
     for j = 1:n1
         
         E(j, j) = b1-2*q1*u(j)*v(j);
         
     end

     E(1, n1)=a1-2*p1*u(n1)*v(n1);
     E(n1, 1)=c1-2*r1*u(1)*v(1);
     %fourth square
     
     for j = 1:n1-1
        
        E(n1+j+1, n1+j) = A1-2*P1*u(j)*v(j);

     end

     for j = 2:n1
         
        E(n1+j-1, n1+j) = C1-2*R1*u(j)*v(j);

     end
     for j = 1:n1
         
         E(n1+j, n1+j) = B1-2*Q1*u(j)*v(j);
         
     end

     E(n1+1, n1+n1)=A1-2*P1*u(n1)*v(n1);
     E(n1+n1, 1+n1)=C1-2*R1*u(1)*v(1);
     
     %second square
     
     for j = 1:n1-1
        
        E(j+1, n1+j) = d1-p1*(u(j)^2+3*v(j)^2);

     end

     for j = 2:n1
         
        E(j-1, n1+j) = f1-r1*(u(j)^2+3*v(j)^2);

     end
     for j = 1:n1
         
         E(j, n1+j) = e1-q1*(u(j)^2+3*v(j)^2);
         
     end

     E(1, n1+n1)=d1-p1*(u(n1)^2+3*v(n1)^2);
     E(n1, 1+n1)=f1-r1*(u(1)^2+3*v(1)^2);
     
     %third square
     for j = 1:n1-1
        
        E(n1+j+1, j) = D1-P1*(v(j)^2+3*u(j)^2);

     end

     for j = 2:n1
         
        E(n1+j-1, j) = F1-R1*(v(j)^2+3*u(j)^2);

     end
     for j = 1:n1
         
         E(n1+j, j) = E1-Q1*(v(j)^2+3*u(j)^2);
         
     end

     E(n1+1, n1)=D1-P1*(v(n1)^2+3*u(n1)^2);
     E(n1+n1, 1)=F1-R1*(v(1)^2+3*u(1)^2);
     
     

     
      eps = inv(E)*f11;
      
      u0 = u+eps(1:n1);
      ep1(:,i) = eps(1:n1);
      ep2(:,i) = eps(n1+1:2*n1);
      v0 = v + eps(n1+1:2*n1);
      cur = u0.^2+v0.^2;
      cur2 = (u0.^2+v0.^2).^2+u01.^2+v01.^2;
      summ = 0;
      summ2 = 0;
      for j = 1:length(cur)
            summ = summ + cur(j);
            summ2 = summ2 + cur2(j);
      end
      integral(i+1) = summ;
      integral2(i+1) = summ2;
      %spectr(i+1) = real(min(eig(E)));
      if i == 10
          anss = E;
      end

      U1(:, i+1) = u0;
      V1(:, i+1) = v0;
      
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
      
%       E1 = diag(ones(n1, 1)*2/3)+diag(ones(n1-1,1)*1/6,1)+diag(ones(n1-1,1)*1/6,-1);
%       E1(1, n1) = 1/6;
%       E1(n1, 1) = 1/6;
%       ans_cur_u = inv(E1)*u01;
%       ans_cur_v = inv(E1)*v01;
      
      cur2 = -u0.*v01+v0.*u01;
      summ1 = 0;
      for j = 1:length(cur)
            summ1 = summ1 + cur2(j);
      end
      integral1(i+1) = summ1;

      

   
    
 end