function u = CoherenceFilter(u,Options)

try
    functionname='CoherenceFilter.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
    addpath([functiondir '/functions2D'])
    addpath([functiondir '/functions3D'])
    addpath([functiondir '/functions'])
catch me
    disp(me.message);
end

% Default parameters
defaultoptions=struct('T',2,'dt',[],'sigma', 1, 'rho', 1, 'TensorType', 1, 'eigenmode',0,'C', 1e-10, 'm',1,'alpha',0.001,'lambda_e',0.02,'lambda_c',0.02,'lambda_h',0.5,'RealDerivatives',false,'Scheme','R','verbose','iter');

if(~exist('Options','var')),
    Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))),
        warning('CoherenceFilter:unknownoption','unknown options found');
    end
end

if(isempty(Options.dt))
    switch lower(Options.Scheme)
      case 'r', Options.dt=0.15;
      case 'o', Options.dt=0.15;
      case 'i', Options.dt=0.15;
      case 's', Options.dt=0.15;
      case 'n', Options.dt=0.15;
      otherwise
        error('CoherenceFilter:unknownoption','unknown scheme');
    end
end
    

% Initialization
dt_max = Options.dt; t = 0;

% In case of 3D use single precision to save memory
if(size(u,3)<4), u=double(u); else u=single(u); end 

% Process time
process_time=tic;

% Show information 
switch lower(Options.verbose(1))
case 'i'
    disp('Diffusion time   Sec. Elapsed');
case 'f'
    disp('Diffusion time   Sec. Elapsed   Image mean    Image variance');
end        


% Anisotropic diffusion main loop
while (t < (Options.T-0.001))
    % Update time, adjust last time step to exactly finish at the wanted
    % diffusion time
    Options.dt = min(dt_max,Options.T-t); t = t + Options.dt;
    tn=toc(process_time);
    switch lower(Options.verbose(1))
    case 'n'
    case 'i'
        s=sprintf('    %5.0f        %5.0f    ',t,round(tn)); disp(s);
    case 'f'
        s=sprintf('    %5.0f        %5.0f      %13.6g    %13.6g ',t,round(tn), mean(u(:)), var(u(:))); disp(s);
    
    end        
   
    if(size(u,3)<4) % Check if 2D or 3D
        % Do a diffusion step
        if(strcmpi(Options.Scheme,'R')&&(Options.eigenmode==0)&&(exist('CoherenceFilterStep2D')==3))
            u=CoherenceFilterStep2D(u,Options);
        else
            u=Anisotropic_step2D(u,Options);
        end
    else
        % Do a diffusion step
        if(strcmpi(Options.Scheme,'R'))
            u=CoherenceFilterStep3D(u,Options);
        else
            u=Anisotropic_step3D(u,Options);
        end
    end
end

function u=Anisotropic_step2D(u,Options)
% Perform tensor-driven diffusion filtering update

% Gaussian smooth the image, for better gradients
usigma=imgaussian(u,Options.sigma,4*Options.sigma);


% Calculate the gradients
switch lower(Options.Scheme)
  case {'r','o','i'}
    ux=derivatives(usigma,'x'); uy=derivatives(usigma,'y');
  case {'s','n'}
    [uy,ux]=gradient(usigma);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end


% Compute the 2D structure tensors J of the image
[Jxx, Jxy, Jyy] = StructureTensor2D(ux,uy,Options.rho);

% Compute the eigenvectors and values of the strucure tensors, v1 and v2, mu1 and mu2
[mu1,mu2,v1x,v1y,v2x,v2y]=EigenVectors2D(Jxx,Jxy,Jyy);

% Gradient magnitude squared
gradA=ux.^2+uy.^2;

% Construct the edge preserving diffusion tensors D = [Dxx,Dxy;Dxy,Dyy]
[Dxx,Dxy,Dyy]=ConstructDiffusionTensor2D(mu1,mu2,v1x,v1y,v2x,v2y,gradA,Options);

% Do the image diffusion
switch lower(Options.Scheme)
  case 'o'
      u=diffusion_scheme_2D_novel(u,Dxx,Dxy,Dyy,Options.dt);
      %u=diffusion_scheme_2D_high_rotation(u,Dxx,Dxy,Dyy,Options.dt,b);
  case 'r'
      u=diffusion_scheme_2D_rotation_invariant(u,Dxx,Dxy,Dyy,Options.dt);
  case 'i'
      u=diffusion_scheme_2D_implicit(u,Dxx,Dxy,Dyy,Options.dt);
  case 's'
      u=diffusion_scheme_2D_standard(u,Dxx,Dxy,Dyy,Options.dt);
  case 'n'
      u=diffusion_scheme_2D_non_negativity(u,Dxx,Dxy,Dyy,Options.dt);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end

function u=Anisotropic_step3D(u,Options)
% Perform tensor-driven diffusion filtering update

% Gaussian smooth the image, for better gradients
usigma=imgaussian(u,Options.sigma,4*Options.sigma);

% Calculate the gradients
ux=derivatives(usigma,'x');
uy=derivatives(usigma,'y');
uz=derivatives(usigma,'z');

% Compute the 3D structure tensors J of the image
[Jxx, Jxy, Jxz, Jyy, Jyz, Jzz] = StructureTensor3D(ux,uy,uz, Options.rho);

% Gradient magnitude squared
gradA=ux.^2+uy.^2+uz.^2;

% Free memory
clear ux; clear uy; clear uz;

% Compute the eigenvectors and eigenvalues of the hessian and directly
% use the equation of Weickert to convert them to diffusion tensors
[Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]=StructureTensor2DiffusionTensor3D(Jxx,Jxy,Jxz,Jyy,Jyz,Jzz,gradA,Options); 

% Free memory
clear J*;

% Do the image diffusion
switch lower(Options.Scheme)
  case 'o'
      u=diffusion_scheme_3D_novel(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
      %u=diffusion_scheme_3D_high_rotation(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 'r'
      u=diffusion_scheme_3D_rotation_invariant(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 'i'
      u=diffusion_scheme_3D_implicit(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 's'
      u=diffusion_scheme_3D_standard(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  case 'n'
      u=diffusion_scheme_3D_non_negativity(u,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,Options.dt);
  otherwise
    error('CoherenceFilter:unknownoption','unknown scheme');
end




        