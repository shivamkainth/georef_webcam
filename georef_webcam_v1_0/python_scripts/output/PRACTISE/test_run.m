vs=1;
os=0; 
cs=1;
rs=0;
rsps=0;
rsms=0;
is=1;
fin_demW='input\6.asc';
fin_folder='output\PRACTISE\';
fin_imfolder='test_run\'; 
if is==1                                                                   %KEEP
    fin_imformat='.jpg';
end
if vs==0 && is==1                                                           %KEEP
    fin_vsformat='.view.asc';
elseif vs==0 && is==0                                                       %KEEP
        fin_viewW='XYZ.view.asc';
end                                                                        %KEEP
if os>0 && is==1                                                            %KEEP
    fin_gcpformat='.gcp.txt';    
end                                                                        %KEEP
cam(:,1)=[76.46626,30.96866]; % altitude will be calculated!
if vs==1                                                                   %KEEP
    buffer_radius=100.0; % in  [m]    
end                                                                        %KEEP
cam(:,2)=[795.8060603386512,-663.6897104589972]; % altitude will be calculated!
cam_off=[209,189.67301929153504]; 
cam_rol=-2; 
cam_foc=0.014;
cam_hei=0.0149; 
cam_wid=0.0223; 
if cs>0                                                                %KEEP
    thres_b_orig=127;                                                      %keep
    movavgwindow=5;
end                                                                %KEEP
if os>1                                                                    %KEEP
    UBD=[50, 50, 50, 3, 100, 100, 100, 0.00010, 0, 0];
    LBD=[-50, -50, -50, -3, -100, -100, -100, -0.00010, 0, 0];
    DDS_R_os=0.2;                                                          %keep
    DDS_MaxEval_os=3000;                                                   %keep
    if os==3                                                               %KEEP
        gcpRMSE_optthres=1;
        gcpRMSE_countthres=10;
    end                                                                    %KEEP
end                                                                       %KEEP
if os>0
    os_MarkSiz=8;
end
cs_MarkSiz=6;
