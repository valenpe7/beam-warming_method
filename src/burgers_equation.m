function burgers_equation
global T t J x un B C
B = 20;
J = 200;
T = 5;
C = 0.8;
t = 0;
nstep = 1;
init;
while t < T
  bwstep;
  plot(x, un, 'x-');
  axis([-B B 0 1]);
  drawnow;
  nstep = nstep + 1;
end

end
  
function init
global dt dx J x u u2 un B C F
dx = 2*B/J;
dt = C*dx/abs(max(u));
x = zeros(J, 1);
for j=1:J
    x(j) = -B + j*dx;
end
F = zeros(J, 1);
u = zeros(J, 1);
u2 = zeros(J, 1);
un = zeros(J, 1);
for j = 1:J
%    if x(j) > 0
%      u(j) = 0.5;
%    else
%      u(j) = -0.5;
%    end
%    if x(j) > 0
     if x(j) < 0
         u(j) = (1 + cos(x(j)/2))/2;
%        u(j) = 1;
     end
%    u(j) = 1/(1+x(j)^2);
end
end

function bwstep
global t dt dx J u un F
t = t + dt;
flux(u, J);
for j = 3:J-2
  if u(j) > 0
    un(j) = u(j) - dt/(2*dx)*(3*F(j) - 4*F(j-1) + F(j-2)) + dt^2/(2*dx^2)*((u(j-1) + u(j))/2*(F(j) - F(j-1)) - (u(j-2) + u(j-1))/2*(F(j-1) - F(j-2)));
  else
    un(j) = u(j) + dt/(2*dx)*(3*F(j) - 4*F(j+1) + F(j+2)) + dt^2/(2*dx^2)*((u(j+1) + u(j+2))/2*(F(j+2) - F(j+1)) - (u(j) + u(j+1))/2*(F(j+1) - F(j)));
  end
end
un(2) = un(3);
un(1) = u(2);
un(J-1) = un(J-2);
un(J) = u(J-1);
u = un;
end

function flux(uu, JJ)
global F k
for k = 1:JJ
    F(k) = uu(k)^2/2;
end
end


