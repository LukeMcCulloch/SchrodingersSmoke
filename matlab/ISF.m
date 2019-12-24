classdef ISF < TorusDEC
% ISF a handle class for simulating incompressible Schroedinger flow.
%
% SYNTAX
%
%   isf = ISF(sizex,sizey,sizez,resx,resy,resz)
%   isf = ISF(sizex,sizey,sizez,res)
%
% DESCRIPTION
%
%   ISF is a subclass of TorusDEC handle class. See TorusDEC for calling
%   constructor.
%
%   To setup, 
%
%     isf = ISF(5,5,5,32,32,32);  % call constructor
%     isf.hbar = 0.1;             % specify Planck constant
%     isf.dt   = 1/24;            % specify time step
%     isf.BuildSchroedinger;      % this command builds coeff for solving
%                                 % Schroedinger equation in Fourier domain
%
%   Useful functions:
%
%     [psi1,psi2] = isf.SchroedingerFlow(psi1,psi2)
%         solves Schroedinger equation for (psi1,psi2) for isf.dt time.
%
%     [psi1,psi2] = isf.Normalize(psi1,psi2)
%         normalizes (psi1,psi2).
%
%     [psi1,psi2] = isf.PressureProject(psi1,psi2)
%         projects (psi1,psi2) to satisfy divergence free constraint.
%
%     [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2,isf.hbar)
%         extracts velocity 1-from from (psi1,psi2)
%
%     [vx,vy,vz] = isf.VelocityOneForm(psi1,psi2)
%         extracts velocity 1-form assuming hbar=1.
%
%
% See also TorusDEC
    
    properties
        hbar             % reduced Planck constant
        dt               % time step
        SchroedingerMask % Fourier coefficient for solving Schroedinger eq
    end
    methods
        
        function obj = ISF(varargin)
        % calls constructor of TorusDEC.
        %
            obj = obj@TorusDEC(varargin{:});
        end
        
        function BuildSchroedinger(obj)
        % builds coefficients in Fourier space.
        %
            nx=obj.resx; ny=obj.resy; nz=obj.resz;
            fac = -4*pi^2*obj.hbar;
            kx = (obj.iix-1-nx/2)/(obj.sizex);
            ky = (obj.iiy-1-ny/2)/(obj.sizey);
            kz = (obj.iiz-1-nz/2)/(obj.sizez);
            lambda = fac*(kx.^2+ky.^2+kz.^2);
            obj.SchroedingerMask = exp(1i*lambda*obj.dt/2);
        end
        
        function [psi1,psi2] = SchroedingerFlow(obj,psi1,psi2)
        % solves Schroedinger equation for dt time.
        %
            psi1 = fftshift(fftn(psi1)); psi2 = fftshift(fftn(psi2));
            psi1 = psi1.*obj.SchroedingerMask;
            psi2 = psi2.*obj.SchroedingerMask;
            psi1 = ifftn(fftshift(psi1)); psi2 = ifftn(fftshift(psi2));
        end
        
        function [psi1,psi2] = PressureProject(obj,psi1,psi2)
        % Pressure projection of 2-component wave function.
        %
            [vx,vy,vz] = obj.VelocityOneForm(psi1,psi2);
            div = obj.Div(vx,vy,vz);
            q = obj.PoissonSolve(div);
            [psi1,psi2] = obj.GaugeTransform(psi1,psi2,-q);
        end
        
        function [vx,vy,vz] = VelocityOneForm(obj,psi1,psi2,hbar)
        % extracts velocity 1-form from (psi1,psi2).
        % If hbar argument is empty, hbar=1 is assumed.
            ixp = mod(obj.ix,obj.resx) + 1;
            iyp = mod(obj.iy,obj.resy) + 1;
            izp = mod(obj.iz,obj.resz) + 1;
            vx = angle(conj(psi1).*psi1(ixp,:,:) ...
                      +conj(psi2).*psi2(ixp,:,:));
            vy = angle(conj(psi1).*psi1(:,iyp,:) ...
                      +conj(psi2).*psi2(:,iyp,:));
            vz = angle(conj(psi1).*psi1(:,:,izp) ...
                      +conj(psi2).*psi2(:,:,izp));
            if nargin<4
                hbar = 1;
            end
            vx = vx*hbar;
            vy = vy*hbar;
            vz = vz*hbar;
        end
        function psi = AddCircle(obj,psi,center,normal,r,d)
        % adds a vortex ring to a 1-component wave function psi.
        % Inputs center, normal, r specify the circle.
        % Input d specify the thickness around the disk to create a boost
        % in phase. Usually d = 5*dx where dx is grid edge length.
            rx = obj.px - center(1);
            ry = obj.py - center(2);
            rz = obj.pz - center(3);
            normal = normal/norm(normal,2);
            alpha = zeros(size(rx));
            z = rx*normal(1) + ry*normal(2) + rz*normal(3);
            inCylinder = rx.^2+ry.^2+rz.^2 - z.^2 < r^2;
            inLayerP = z> 0 & z<= d/2 & inCylinder;
            inLayerM = z<=0 & z>=-d/2 & inCylinder;
            alpha(inLayerP) = -pi*(2*z(inLayerP)/d - 1);
            alpha(inLayerM) = -pi*(2*z(inLayerM)/d + 1);
            psi = psi.*exp(1i*alpha);
        end
    end
    methods (Static)
        function [psi1,psi2] = GaugeTransform(psi1,psi2,q)
        % multiplies exp(i*q) to (psi1,psi2)
        %
            eiq = exp(1i*q);
            psi1 = psi1.*eiq;
            psi2 = psi2.*eiq;
        end
        function [sx,sy,sz] = Hopf(psi1,psi2)
        % extracts Clebsch variable s=(sx,sy,sz) from (psi1,psi2)
        %
            a = real(psi1);
            b = imag(psi1);
            c = real(psi2);
            d = imag(psi2);
            sx = 2*(a.*c + b.*d);
            sy = 2*(a.*d - b.*c);
            sz = a.^2 + b.^2 - c.^2 - d.^2;
        end
        function [psi1,psi2] = Normalize(psi1,psi2)
        % normalizes (psi1,psi2).
        %
            psi_norm = sqrt(abs(psi1).^2 + abs(psi2).^2);
            psi1 = psi1./psi_norm;
            psi2 = psi2./psi_norm;
        end
    end
end