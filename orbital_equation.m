function [pos] = orbital_equation(m0,m1,G,mu, r0, v0)

d0 = norm(r0); %Initia distance between the two bosies

%Orbital time period
tp = 2 * pi * ( ...
    sqrt(d0^3 / mu) ...
);

s0 = norm(v0); %initial speed


%Initial angular momentum vector
h0 = cross(r0,v0);

%Calculate eccentricity
e = sqrt(1 +  ...
    (norm(h0)^2 / mu^2) * (s0^2 - (2*mu/d0))...
    );

%Calculate eccentricity vector
eVec = (1/mu) * ( ...
    r0*(s0^2 - mu/d0) - v0*(dot(r0,v0)));


% Calculate the Mechanical energy
E = (s0^2 / 2) - (mu / d0);

%Calculate the semi-major axis
a = - ( ...
    mu / (2*E) ...
    );

%Calculate the inclination:
i = acos( ...
    h0(3) / norm(h0));

%Define a variable \hat{k} that is equivalent to 0i + 0j + 1k
Khat = [0, 0, 1];

%Calculate the node line unitary vector 
N = (cross(Khat, h0))/(norm(cross(Khat, h0)));

%Calculate the RAAN
if N(2) < 0
    omega = (2*pi) - acos(N(1) / norm(N));
else
    omega = acos(N(1) / norm(N));
end

%Calculate the argument of periapsis
if dot(h0, cross(N, eVec))<0
     w = 2*pi - acos(dot(N, eVec)/(norm(N)*e));
else
    w = acos(dot(N, eVec)/(norm(N)*e));
end

%Calculate the true anomaly:
if dot(h0, (cross(eVec, r0))) <0
    theta = 2*pi - acos((dot(r0,eVec))/d0*e);
else
    theta = acos((dot(r0,eVec))/d0*e);
end

%Calculate semi-parameter
p = a*(1-e^2);

%Calculate the orbital equation:
oe = p / (1+e*cos(theta));

%% Plotting:

thetas = linspace(0,2*pi,100);

position = [(p*cos(thetas))./(1+e*cos(thetas)); (p*sin(thetas))./(1+e*cos(thetas)); zeros(length(thetas),1).' ].';

rdot = sqrt( ...
    (mu/p)*e*sin(thetas) ... 
    ).';

velocity = [-sqrt(mu/p)*sin(thetas);sqrt(mu/p)*(e+cos(thetas)) ; zeros(length(thetas),1).' ].';

rotMatrix = [cos(omega)*cos(w) - sin(omega)*sin(w)*cos(i), -cos(omega)*sin(w) - sin(omega)*cos(w)*cos(i), sin(omega)*sin(i);
       sin(omega)*cos(w) + cos(omega)*sin(w)*cos(i), -sin(omega)*sin(w) + cos(omega)*cos(w)*cos(i), -cos(omega)*sin(i);
       sin(omega)*sin(i), cos(omega)*sin(i), cos(i)];

pos = real(position*rotMatrix);

end

