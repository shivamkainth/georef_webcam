% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Karsten Schulz, Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xopt, Fopt] = DDS_PRACTISE(FUN, X, LB, UB, R, MaxEval, Var1, ...
    Var2, Var3, Var4, Var5)
%   Name:       DDS_PRACTISE
%   Purpose:    Dynamically Dimensioned Search (DDS) Algorithms searching 
%               for a global minimum of any objective function (Tolson & 
%               Shoemaker, 2007) adapted to PRACTISE
%
%   Output:     Fopt (Value of the objective function in the optimum)
%               Xopt (Parameters for optimal objective function)
%   Input:      FUN (Handle to the objective Function)
%               X (Initial guess for the optimal parameters)
%               LB (Lower boundary of the parameter space)
%               UB (Upper boundary of the parameter space)
%               R (Neighborhood pertubation size)
%               MaxEval (Maximum number of function evaluations)
%           for GCPs
%               Var1=dem (DEM raster from ASCII-file [cf. ESRI ArcGIS])
%               Var2=header (values of DEM ASCII-file header [cf. ESRI ArcGIS]:
%                 ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
%               Var3=gcpW
%               Var4=pix_c
%               Var5=pix_r
%           for NDSI
%               Var1=classified photo map
%               Var2=NDSI values of satellite image (same size as Var1)
%               Var3=name of the objective function (MBCT, BCT1, BCT2)
%               Var4=probabilities of photograph classification (de)activated
%%%   --optional  Var5=number of digits (Xopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: DDS parameters
%   specified in parameter file
%
% Step 2: initial evaluation of the function
if size(Var1,1)*size(Var1,2)~=size(Var2,1)*size(Var2,2)
    [F, gcpP] = feval(FUN, X, Var1, Var2, Var3, Var4, Var5); 
%     % plot initial projected GCP %uncomment if interested in optimisation pathway
%     plot(gcpP(1,:)+0.5, gcpP(2,:)+0.5, 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerSize',8)
%     drawnow
elseif size(Var1,1)*size(Var1,2)==size(Var2,1)*size(Var2,2)
    [F] = feval(FUN, X, Var1, Var2, Var3, Var4); 
else
    error('Wrong number of input arguments in function DDS_PRACTISE.')
end
disp('The initial objective function value is')
disp(['F0 = ', num2str(F)])
%
count   = 1;
disp('The number of DDS runs evaluated is')
disp(count)
Fopt    = F;
Xopt    = X;
quit    = 0;
%
while(quit==0)
%    
    % Step 3: Randomly select J of the N variables
    Pi = 1-log(count)/log(MaxEval);
    Pn = rand(size(X));
    J=find(Pn<Pi);
    if (length(J)==0)
        J=floor(rand(1)*length(X)+1);
    end
    %
    % Step 4: calculate updated xnew
    sig    = R*(UB-LB).*randn(size(X));
    Xnew   = Xopt;
    Xnew(J)= Xnew(J)+sig(J);
    %
    i      = find(Xnew<LB);
    Xnew(i)= LB(i)+(LB(i)-Xnew(i));
    i1     = find(Xnew(i)>UB(i));
    Xnew(i1)=LB(i1);
    %
    z       = find(Xnew>UB);
    Xnew(z) = UB(z) - (Xnew(z)-UB(z));
    z1      = find(Xnew(z)<LB(z));
    Xnew(z1)= UB(z1);
    %
    % Step 5: 
    count   = count+1;
    if size(Var1,1)*size(Var1,2)~=size(Var2,1)*size(Var2,2)
        %%
%         disp(['*', num2str(count), '*'])
%         sprintf('%16.7f',Xnew)
        %%
        [F,gcpP] = feval(FUN, Xnew, Var1, Var2, Var3, Var4, Var5);
    elseif size(Var1,1)*size(Var1,2)==size(Var2,1)*size(Var2,2)
%         if exist('Var5','var')
%             Xnew = round(10^Var5*Xnew)/10^Var5;
%         end
        [F] = feval(FUN, Xnew, Var1, Var2, Var3, Var4); 
    end
    if (F<Fopt)
            Fopt    = F;
            disp(['Fopt = ', num2str(F)])
            Xopt    = Xnew;
%         if size(Var1,1)*size(Var1,2)~=size(Var2,1)*size(Var2,2) %uncomment if interested in optimisation pathway
%             % plot improved projected GCPs
%             plot(gcpP(1,:)+0.5, gcpP(2,:)+0.5, '.', 'MarkerEdgeColor', [1 0 0], 'MarkerSize',8)
%             drawnow
%         end
    end
    %
    % Step 6: Stop Criteria
    if size(Var1,1)*size(Var1,2)~=size(Var2,1)*size(Var2,2) && ...
     (mod(count,500)==0)
        disp(count)
    elseif size(Var1,1)*size(Var1,2)==size(Var2,1)*size(Var2,2) && ...
     (mod(count,50)==0)
        disp(count)
    end    
    if(count==MaxEval)
        quit = 1;
    end    
end