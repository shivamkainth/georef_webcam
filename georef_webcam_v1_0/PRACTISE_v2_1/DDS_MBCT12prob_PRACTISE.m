% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   11/2013
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ObjVal] = DDS_MBCT12prob_PRACTISE(X, Var1, Var2, ObjFunName, Prob)
%   Name:       DDS_MBCT_PRACTISE
%   Purpose:    DDS calls MBCT as objective function to calculate "1-" an 
%               addition of the binary contingency  tables (ABCT) of  
%               observed and modelled snow cover and non snow covered 
%               grid elements
%   
%   Output:     ObjVal ("1-" multiplication of performance measures which  
%               are in optimum=1)
%   Input:      X=estimated NDSI threshold
%               Var1=classified photo map
%               Var2=NDSI values of satellite image (same size as Var1)
%               ObjFunName=name of the objective function (MBCT, BCT1, BCT2)
%               Prob=probabilities of photograph classification (de)activated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create snow map using the satellite image and the NDSI threshold
satBBmapOpt=Var2;
satBBmapOpt(satBBmapOpt>=X)=1; %snow (including threshold value)
satBBmapOpt(satBBmapOpt<X)=0; %no snow

% a(os+ms, hits) b(ons+ms, fails) c(os+mns, misses) d(ons+mns, zero)
if Prob==0
    a=sum(Var1+satBBmapOpt==2);
    b=sum(Var1-satBBmapOpt==-1);
    c=sum(Var1-satBBmapOpt==1);
    d=sum(Var1+satBBmapOpt==0);
elseif Prob==1
    i=find(round(Var1)+satBBmapOpt==2);
    a=sum(Var1(i));
    i=find(round(Var1)-satBBmapOpt==-1);
    b=sum(1-Var1(i));
    i=find(round(Var1)-satBBmapOpt==1);
    c=sum(Var1(i));
    i=find(round(Var1)+satBBmapOpt==0);
    d=sum(1-Var1(i));
    clear i
end
n=a+b+c+d;

F1=(a+d)/n;
F2=a/(a+b+c);
% minimisation problem
if strcmp(ObjFunName, 'MBCT') 
    ObjVal=(1-F1)*(1-F2);
elseif strcmp(ObjFunName, 'BCT1') 
    ObjVal=(1-F1);
elseif strcmp(ObjFunName, 'BCT2') 
    ObjVal=(1-F2);
end