function [g,h] = constraints(L)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

E = [0 -1;
     1  0];
base = 174;
pb = (sqrt(3)/2)*base;
r = 50;
platform = 15;
L1 = L(1);
L2 = L(2);
pp = (sqrt(3)/2)*platform;
A = [-(base/2), -(pb/3);
       base/2, -(pb/3);
       0, (2*pb)/3];
i = 1;
for alpha = 0: pi/4: (7*pi)/4

    
    Boundary(i,:) = [r*cos(alpha), r*sin(alpha)];
    i = i+1;
    
end

phi = pi/6;
for j = 1: 8
    c11(j, :) = [Boundary(j,1)+pp*cos(phi), Boundary(j,2)+pp*sin(phi)];
    c12(j, :) = [Boundary(j,1), Boundary(j,2)-pp];
    c13(j, :) = [Boundary(j,1)-pp*cos(phi), Boundary(j,2)+pp*sin(phi)];
        
    c21(j, :) = [Boundary(j,1)+pp*cos(phi), Boundary(j,2)-pp*sin(phi)];
    c22(j, :) = [Boundary(j,1)-pp*cos(phi), Boundary(j,2)-pp*sin(phi)];
    c23(j, :) = [Boundary(j,1), Boundary(j,2)+pp];
    
    
end
m = 1;
for k = 1: 8
    A1C11 = [sqrt(((A(1,1)-c11(k,1))^2)+((A(1,2)-c11(k,2))^2))];
    A1C21 = [sqrt(((A(1,1)-c21(k,1))^2)+((A(1,2)-c21(k,2))^2))];
    
    A2C12 = [sqrt(((A(2,1)-c12(k,1))^2)+((A(2,2)-c12(k,2))^2))];
    A2C22 = [sqrt(((A(2,1)-c22(k,1))^2)+((A(2,2)-c22(k,2))^2))];
    
    A3C13 = [sqrt(((A(3,1)-c13(k,1))^2)+((A(3,2)-c13(k,2))^2))];
    A3C23 = [sqrt(((A(3,1)-c23(k,1))^2)+((A(3,2)-c23(k,2))^2))];    
    
    B1y = ((A1C11)^2 + (L1)^2 - (L2)^2)/(2*A1C11);
    B1x = sqrt((L1)^2 - (B1y)^2);
    
    
    B1(m, :) = [B1x, B1y];
    B1(m+1, :) = [B1x, B1y];
    
    B2y = ((A2C12)^2 + (L1)^2 - (L2)^2)/(2*A2C12);
    B2x = sqrt((L1)^2 - (B2y)^2);
   
    
    B2(m, :) =  [B2x, B2y];
    B2(m+1, :) = [B2x, B2y];
    
    B3y = ((A3C13)^2 + (L1)^2 - (L2)^2)/(2*A3C13);
    B3x = sqrt((L1)^2 - (B3y)^2);
    
    
    B3(m, :) =  [B3x, B3y];
    B3(m+1, :) = [B3x, B3y];
    m = m+2;
end
    
for j = 1 : 16
    A1B1_vector = [A(1,1)-B1(j,1), A(1,2)-B1(j,2)];
    A1B1 = sqrt((A(1,1)-B1(j,1))^2+(A(1,2)-B1(j,2))^2);
    u1(j, :) = [A1B1_vector/A1B1];
    
    A2B2_vector = [A(2,1)-B2(j,1), A(2,2)-B2(j,2)];
    A2B2 = sqrt((A(2,1)-B2(j,1))^2+(A(2,2)-B2(j,2))^2);
    u2(j, :) = [A2B2_vector/A2B2];
    
    A3B3_vector = [A(3,1)-B3(j,1), A(3,2)-B3(j,2)];
    A3B3 = sqrt((A(3,1)-B3(j,1))^2+(A(3,2)-B3(j,2))^2);
    u3(j, :) = [A3B3_vector/A3B3];
end
p = 1;
for j = 1 : 8
    B1C1(p) = [sqrt(((B1(p,1)-c11(j,1))^2)+((B1(p,2)-c11(j,2))^2))];
    B1C1(p+1) = [sqrt(((B1(p+1,1)-c21(j,1))^2)+((B1(p+1,2)-c21(j,2))^2))];
    B1C1_vector(p, :) = [B1(p,1)-c11(j,1), B1(p,2)-c11(j,2)];
    B1C1_vector(p+1, :) = [B1(p+1,1)-c21(j,1), B1(p+1,2)-c21(j,2)];
    
                   
    B2C2(p) = [sqrt(((B2(p,1)-c12(j,1))^2)+((B2(p,2)-c12(j,2))^2))];
    B2C2(p+1) = [sqrt(((B2(p+1,1)-c22(j,1))^2)+((B2(p+1,2)-c22(j,2))^2))];
    B2C2_vector(p, :) = [B2(p,1)-c11(j,1), B2(p,2)-c11(j,2)];
    B2C2_vector(p+1, :) = [B2(p+1,1)-c21(j,1), B2(p+1,2)-c21(j,2)];
    
                           
    B3C3(p) = [sqrt(((B3(p,1)-c13(j,1))^2)+((B3(p,2)-c13(j,2))^2))];
    B3C3(p+1) = [sqrt(((B3(p+1,1)-c23(j,1))^2)+((B3(p+1,2)-c23(j,2))^2))];
    B3C3_vector(p, :) = [B3(p,1)-c11(j,1), B3(p,2)-c11(j,2)];
    B3C3_vector(p+1, :) = [B3(p+1,1)-c21(j,1), B3(p+1,2)-c21(j,2)];
    p = p+2;
end
for j = 1 : 16
    v1(j, :) = B1C1_vector(j, :)/B1C1(j);
    v2(j, :) = B2C2_vector(j, :)/B2C2(j);
    v3(j, :) = B3C3_vector(j, :)/B3C3(j);
end
q = 1;
 for j = 1 : 8
       C1P_vector(q, :) = [c11(j,1)-Boundary(j,1), c11(j,2)-Boundary(j,2)];
       C1P_vector(q+1, :) = [c21(j,1)-Boundary(j,1), c21(j,2)-Boundary(j,2)];
       
       C2P_vector(q, :) = [c12(j,1)-Boundary(j,1), c12(j,2)-Boundary(j,2)];
       C2P_vector(q+1, :) = [c22(j,1)-Boundary(j,1), c22(j,2)-Boundary(j,2)];
       
       C3P_vector(q, :) = [c13(j,1)-Boundary(j,1), c13(j,2)-Boundary(j,2)];
       C3P_vector(q+1, :) = [c23(j,1)-Boundary(j,1), c23(j,2)-Boundary(j,2)];
       
       C1P(q) = [sqrt(((c11(j,1)-Boundary(j,1))^2)+((c11(j,2)-Boundary(j,2))^2))];
       C1P(q+1) = [sqrt(((c21(j,1)-Boundary(j,1))^2)+((c21(j,2)-Boundary(j,2))^2))];
       
       C2P(q) = [sqrt(((c12(j,1)-Boundary(j,1))^2)+((c12(j,2)-Boundary(j,2))^2))];
       C2P(q+1) = [sqrt(((c22(j,1)-Boundary(j,1))^2)+((c22(j,2)-Boundary(j,2))^2))];
       
       C3P(q) = [sqrt(((c13(j,1)-Boundary(j,1))^2)+((c13(j,2)-Boundary(j,2))^2))];
       C3P(q+1) = [sqrt(((c23(j,1)-Boundary(j,1))^2)+((c23(j,2)-Boundary(j,2))^2))];
       q = q+2;
 end
 for j = 1: 16
    k1(j, :) = C1P_vector(j, :)/C1P(j);
    k2(j, :) = C2P_vector(j, :)/C2P(j);
    k3(j, :) = C3P_vector(j, :)/C3P(j);
 end
 
%  s = 1;
for j = 1 : 16 
    J1 = [-pp*(v1(j, :))*E*(k1(j, :))', v1(j, :);
             -pp*(v2(j, :))*E*(k2(j, :))', v2(j, :);
             -pp*(v3(j, :))*E*(k3(j, :))', v3(j, :)];
    J2 = [L1*v1(j, :)*E*(u1(j, :))', 0, 0;
             0, L1*v2(j, :)*E*(u2(j, :))', 0;
             0, 0, L1*v3(j, :)*E*(u3(j, :))'];
    J = J1\J2;
    c_num = cond(J);
    inv_condition(j) = inv(c_num);
end

    g1 = 0.1-inv_condition;
    g2 = 150-L1+L2;
    
    g = [g1, g2];
    h = [];
end

