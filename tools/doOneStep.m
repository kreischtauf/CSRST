function [Ex Ey Hz] = doOneStep(Lx, Ly, Ex, Ey, Hz, dt, dx, dy)

[Ex Ey Hz] = firstStep(Lx, Ex, Ey, Hz, dt, dx, dy);
[Ex Ey Hz] = secondStep(Ly, Ex, Ey, Hz, dt, dx, dy);