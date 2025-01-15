function [Jxx, Jxy, Jyy]=StructureTensor2D(ux,uy,rho)

% J(grad u_sigma)
Jxx = sum(ux.^2,3)/size(ux,3);
Jxy = sum(ux.*uy,3)/size(ux,3);
Jyy = sum(uy.^2,3)/size(ux,3);

% Do the gaussian smoothing of the structure tensor
Jxx = imgaussian(Jxx,rho,6*rho);
Jxy = imgaussian(Jxy,rho,6*rho);
Jyy = imgaussian(Jyy,rho,6*rho);