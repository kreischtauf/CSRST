function [Ex Ey Hz] = secondStep(Ly, Ex, Ey, Hz, dt, dx, dy)

eps = 8.854187817e-12;
mu = 1.256637061e-06;
csqr = 1 / (mu * eps);
[Ny Nx] = size(Hz);

dtepsdx = dt / (2 * eps * dx);
dtepsdy = dt / (2 * eps * dy);
dtmudx = dt / (2 * mu * dx);
dtmudy = dt / (2 * mu * dy);

tmp = Ex;
tmp(2:Ny,:) = tmp(2:Ny,:) + dtepsdy * (Hz(2:Ny,:) - Hz(1:Ny-1,:)) - ...
    dtmudx * dtepsdy * (Ey(2:Ny,2:Nx+1) - Ey(2:Ny,1:Nx) - ...
                        Ey(1:Ny-1,2:Nx+1) + Ey(1:Ny-1,1:Nx));

Ex = Ly \ tmp;

tmp = Ey;
Ey(:,2:Nx) = Ey(:,2:Nx) - dtepsdx * (Hz(:,2:Nx) - Hz(:,1:Nx-1));
Hz = Hz + dtmudy * (Ex(2:Ny+1,:) - Ex(1:Ny,:)) - ...
     dtmudx * (tmp(:,2:Nx+1) - tmp(:,1:Nx));