function R = calculateRotationalAffine(alpha, beta, gamma)

Ralpha = [1 1 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Rbeta = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rgamma = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

R = Ralpha*Rbeta*Rgamma;

end