function advection_equation
global T J B C a
B = 20;
J = 200;
T = 50;
C = 1;
a = 1;
solve;
end

function solve
global T t dt x u un a B
t = 0;
init;
while t < T
  if a > 0
      beam_warming1;
  else
      beam_warming2;
  end
  plot(x, un, 'x-');
  axis([-B B 0 1]);
  drawnow;
  u = un;
  t = t + dt;
end
end

function init
global dt dx J x u un B C a
dx = 2*B/J;
dt = C*dx/abs(a);
x = zeros(J, 1);
for j = 1:J
    x(j) = -B + j*dx;
end
u = zeros(J, 1);
un = zeros(J, 1);
for j = 1:J
    u(j) = exp(-(x(j)/4)^2);
%    if abs(x(j)) < pi
%       u(j) = (1 + cos(x(j)))/2;
%       u(j) = 1;
%    end
end
end

function beam_warming1
global dt dx J u un a
for j = 3:J
    un(j) = u(j) - a*dt*(3*u(j) - 4*u(j-1) + u(j-2))/(2*dx) + a^2*dt^2*(u(j) - 2*u(j-1) + u(j-2))/(2*dx^2);
end
un(1) = un(J-1);
un(2) = un(J);
end

function beam_warming2
global dt dx J u un a
for j = 1:J-2
    un(j) = u(j) + a*dt*(3*u(j) - 4*u(j+1) + u(j+2))/(2*dx) + a^2*dt^2*(u(j) - 2*u(j+1) + u(j+2))/(2*dx^2);
end
un(J-1) = un(1);
un(J) = un(2);
end
