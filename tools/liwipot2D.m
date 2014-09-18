function [Ex, Ey] = liwipot2D(Ekin, dx, dy, dz, xrange, yrange, zrange)

gamma = Ekin / 0.511 + 1;
beta = [sqrt(1.0 - 1.0 / gamma**2);
        0;
        0];

fh = fopen('liwipot2.dat', 'w');
Ex = zeros(size(yrange,2), size(xrange,2));
Ey = zeros(size(yrange,2), size(xrange,2));
for j = 1:size(yrange,2)
  ry = dy * yrange(j);
  for i = 1:size(xrange,2)
    rx = dx * xrange(i);

    rz = dz * zrange;
    Nz = size(zrange,2);
    one = ones(1,Nz);
    cdt = gamma**2 * (rx * beta(1) + ...
                      sqrt(rx**2 + (ry**2 + rz.**2) / gamma**2));
    rho = -beta * cdt;
    n = [rx*one;ry*one;rz] - rho;
    n = n ./ (ones(3,1) * sqrt(dot(n,n)));

    E = (n - beta*one) ./ (ones(3,1) * (gamma**2 * (1 - beta'*n).^3 .* cdt.**2));
    Ex(j, i) = sum(E(1,:)) / Nz;
    Ey(j, i) = sum(E(2,:)) / Nz;

    [v, k] = min(abs(zrange));
    fprintf(fh, '%15.10f %15.10f %15.10f %15.10f %15.10f\n', ...
            rx, ry, Ex(j,i), Ey(j,i), cdt(k));
  end
end

fclose(fh);