% get a comprehensive Def- Pot. and write scattering matrixes in .mat file

clearvars

save_filename = 'scattering_parameters_HH24.mat' ;

EPC_filename = 'HH24_DeformationPotentials_shifted_data' ; % 'C8BTBTC8_triclinic_DeformationPotentials_data.mat';% 'ROQSAT_G100_DeformationPotentials_data_modified.mat' ; % 'ROQSAT_G100_DeformationPotentials_data.mat';%'rubrene_ortho_02tag_443mesg_DeformationPotentials_data.mat' ; % 'phBTBT10_ROQSAT_DeformationPotentials_data_1mode_removed.mat' ; % phBTBT10_herringbone_DeformationPotentials_4modes_removed.mat' ; %  'phBTBT10_ROQSAT_DeformationPotentials_data.mat'; %'Anthra_1_DeformationPotentials_data.mat';

phononDOS_filename = 'total_dos_HH24.dat' ; % HH3_total_dos 'total_dos_HH24.dat' ; %  'C8BTBTC8_P-1_total_dos.dat'; % 'phBTBT10_herringbone_total_dos.dat' ; % 'phBTBT10_ROQSAT_VdW_total_dos.dat';


f_cut = 8 ; %7.7; % HH24 5.2

rho_mass_density = 7.4e3 ; % 7.7e3; % Ni3In 9.27e3; HH24 7.4e3 7.7e3  HH3 

ADP = 'yes' ; %'yes' or 'no'
us_sound_vel = 2e3;

% optional, if provided the code shrinks the matrixes to fit the number of bands; particularly usefule for metals with narrow bands
% they contain the indexes of the bands of interest in the DFT band structure
bands_transp = 90:96; % 89:98 ; % 90:97; % 10:12; HH24 for EF -0.5:0.5 bands are 88-100
% bands_transp_CB
% bands_transp_VB
%


load(EPC_filename)

bands_type_CB =  1:size(DP_bw(1,1).CB,1);

n_CB = max(size(DP_bw(1,1).CB)) ;

n_m = size(DP_bw,1);
n_q = size(DP_bw,2);

n_process_odp_CB = n_m; % these should be n_m + length of n_m_intra_array
n_process_ivs_CB = n_m;


if strcmp(semiconductor,'yes') % the variable semiconductor is saved in the DefPot .mat file
    bands_type_VB =  1:size(DP_bw(1,1).VB,1);
    n_VB = max(size(DP_bw(1,1).VB)) ;
    
    n_process_odp_VB = n_m;
    n_process_ivs_VB = n_m;
else
    semiconductor = 'no';
end



imported = importdata(phononDOS_filename);
DOS_ph = imported.data;



fr = DOS_ph(DOS_ph(:,1)<f_cut,1);
dosr = DOS_ph(DOS_ph(:,1)<f_cut,2);


D_odp_e_M = zeros(n_process_odp_CB,n_CB);
hbar_w_odp_e_M = zeros(n_process_odp_CB,n_CB);

D_ivs_e_M = zeros(n_CB,n_CB,n_process_ivs_CB);
hbar_w_ivs_e_M = zeros(n_CB,n_CB,n_process_ivs_CB);

if strcmp(semiconductor,'yes')
    D_odp_h_M = zeros(n_process_odp_VB,n_VB);
    hbar_w_odp_h_M = zeros(n_process_odp_VB,n_VB);
    
    D_ivs_h_M = zeros(n_VB,n_VB,n_process_ivs_VB);
    hbar_w_ivs_h_M = zeros(n_VB,n_VB,n_process_ivs_VB);
end


for i_b = 1:n_CB
    
    for im = n_m :-1:1    
        for iq = n_q :-1:1

        f_t = DP_bw(im,iq).freq;

        DOS_t= interp1(fr,dosr,f_t, 'makima');
        if DOS_t < 0 || f_t < 0
            DOS_t = 0;
        end
        DOS_M(im,iq) = DOS_t;

       DP_temp(im,iq) = DP_bw(im,iq).CB(i_b,i_b) ;
        
        end
    end


    for im = n_m:-1:1
        for iq = n_q:-1:1
            freq_temp(im,iq) = DP_bw(im,iq).freq ;
        end
    end
    w_v2_array = ( (sum(freq_temp.*DOS_M,2))./sum(DOS_M,2) )*4.136e-3 ; % in eV



     % D_odp_e_M(1,i_b) = ...
     %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

     D_odp_e_M(:,i_b) = ...
         sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;

     hbar_w_odp_e_M(:,i_b) = w_v2_array;
            
     for i_b2 = 1:n_CB
        if ne(i_b,i_b2)
            for im = n_m:-1:1    
                for iq = n_q:-1:1
                    DP_temp(im,iq) =  DP_bw(im,iq).CB(i_b,i_b2) ;
                end
            end

            D_ivs_e_M(i_b,i_b2,:) = ...
                sqrt(...
         sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
         ) *1e10;
            % D_ivs_e_M(i_b,i_b2,1) = ...
            %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;

            hbar_w_ivs_e_M(i_b,i_b2,:) = w_v2_array;
        end
    end

end

if strcmp(semiconductor,'yes')
    for i_b = 1:n_VB
        
        for im = n_m:-1:1    
            for iq = n_q:-1:1
    
            DP_temp(im,iq) = DP_bw(im,iq).VB(i_b,i_b) ;
            
            end
        end
    
         D_odp_h_M(:,i_b) = ...
             sqrt(...
             sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
             ) *1e10;
             % D_odp_h_M(1,i_b) = ...
             %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;
    
         hbar_w_odp_h_M(:,i_b) = w_v2_array;
                
         for i_b2 = 1:n_VB
            if ne(i_b,i_b2)
                for im = n_m:-1:1    
                    for iq = n_q:-1:1
                        DP_temp(im,iq) =  DP_bw(im,iq).VB(i_b,i_b2) ;
                    end
                end
    
                D_ivs_h_M(i_b,i_b2,:) = ...
                    sqrt(...
             sum(DP_temp.^2.*DOS_M,2)./(sum(DOS_M,2))...
             ) *1e10;
                % D_ivs_h_M(i_b,i_b2,1) = ...
                %     sqrt(sum(sum(DP_temp.^2.*DOS_M))/(sum(sum(DOS_M)))) *1e10;
    
                hbar_w_ivs_h_M(i_b,i_b2,:) = w_v2_array;
            end
        end
   end
end

if min(hbar_w_odp_e_M)<0
    hw_f = min(hbar_w_odp_e_M(hbar_w_odp_e_M>0)) / 2;
    hbar_w_odp_e_M(hbar_w_odp_e_M<0) = hw_f;
    hw_f = min(hbar_w_ivs_e_M(hbar_w_ivs_e_M>0)) / 2;
    hbar_w_ivs_e_M(hbar_w_ivs_e_M<0) = hw_f;
    if strcmp(semiconductor,'yes')
          hw_f = min(hbar_w_odp_h_M(hbar_w_odp_h_M>0)) / 2;
          hbar_w_odp_h_M(hbar_w_odp_h_M<0) = hw_f;
          hw_f = min(hbar_w_ivs_h_M(hbar_w_ivs_h_M>0)) / 2;
          hbar_w_ivs_h_M(hbar_w_ivs_h_M<0) = hw_f;
    end      
end



%
k_s = 100;
k_inf = 100; 
Z_i = 1;
%
%
%
%
%
if strcmp(ADP,'yes')
    if strcmp(semiconductor,'yes')
        if exist('bands_transp_CB','var') && isempty(bands_transp_CB) == 0
            if min(CB_idx)>min(bands_transp_CB) || max(CB_idx)<max(bands_transp_CB)
                disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extracted for a larger energy window')
            else
                 [v,p] = intersect(CB_idx,bands_transp_CB);
                 temp = D_adp_CB(p,p) ;
                 D_adp_CB = temp;
            end    
        end
        if exist('bands_transp_VB','var') && isempty(bands_transp_VB) == 0
            if min(VB_idx)>min(bands_transp_VB) || max(VB_idx)<max(bands_transp_VB)
                disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extracted for a larger energy window')
            else
                 [v,p] = intersect(VB_idx,bands_transp_VB);
                 temp = D_adp_VB(p,p) ;
                 D_adp_VB = temp;
            end
        end
        D_adp_h_M = D_adp_VB;
    else
        if exist('bands_transp','var') && isempty(bands_transp) == 0
            if min(CB_idx)>min(bands_transp) || max(CB_idx)<max(bands_transp)
                disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extraacted for a larger energy window')
            else
                 [v,p] = intersect(CB_idx,bands_transp);
                 temp = D_adp_CB(p,p) ;
                 D_adp_CB = temp;
            end    
        end
    end
    D_adp_e_M = D_adp_CB;
end
%
%
%
NNan = isnan(D_odp_e_M);
D_odp_e_M(NNan) = 0;
NNan = isnan(hbar_w_odp_e_M);
hbar_w_odp_e_M(NNan) = 0;
NNan = isnan(D_ivs_e_M);
D_ivs_e_M(NNan) = 0;
NNan = isnan(hbar_w_ivs_e_M);
hbar_w_ivs_e_M(NNan) = 0;
if strcmp(semiconductor,'yes')
    NNan = isnan(D_odp_h_M);
    D_odp_h_M(NNan) = 0;
    NNan = isnan(hbar_w_odp_h_M);
    hbar_w_odp_h_M(NNan) = 0;
    NNan = isnan(D_ivs_h_M);
    D_ivs_h_M(NNan) = 0;
    NNan = isnan(hbar_w_ivs_h_M);
    hbar_w_ivs_h_M(NNan) = 0;
end

%
%
%
if strcmp(semiconductor,'yes')
    if exist('bands_transp_CB','var') && isempty(bands_transp_CB) == 0
        if min(CB_idx)>min(bands_transp_CB) || max(CB_idx)<max(bands_transp_CB)
            disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extracted for a larger energy window')
        else
             [v,p] = intersect(CB_idx,bands_transp_CB);
             temp = D_ivs_e_M(p,p,:) ;
             D_ivs_e_M = temp;
             temp = D_odp_e_M(:,p) ;
             D_odp_e_M = temp;
             temp = hbar_w_ivs_e_M(p,p,:) ;
             hbar_w_ivs_e_M = temp;
             temp = hbar_w_odp_e_M(:,p) ;
             hbar_w_odp_e_M = temp;
        end    
    end
    if exist('bands_transp_VB','var') && isempty(bands_transp_VB) == 0
        if min(VB_idx)>min(bands_transp_VB) || max(VB_idx)<max(bands_transp_VB)
            disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extracted for a larger energy window')
        else
             [v,p] = intersect(VB_idx,bands_transp_VB);
             temp = D_ivs_h_M(p,p,:) ;
             D_ivs_h_M = temp;
             temp = D_odp_h_M(:,p) ;
             D_odp_h_M = temp;
             temp = hbar_w_ivs_h_M(p,p,:) ;
             hbar_w_ivs_h_M = temp;
             temp = hbar_w_odp_h_M(:,p) ;
             hbar_w_odp_h_M = temp;
        end    
    end
    if strcmp(ADP,'yes')
        save(save_filename, "Z_i", 'k_inf', 'k_s', 'hbar_w_ivs_h_M', 'D_ivs_h_M'...
        , "hbar_w_odp_h_M", "D_odp_h_M", "hbar_w_ivs_e_M", "hbar_w_odp_e_M" ...
        , "D_adp_e_M", "D_adp_h_M", "us_sound_vel"...
        , "D_odp_e_M", "D_ivs_e_M", "rho_mass_density", "bands_type_VB", "bands_type_CB")
    else
        save(save_filename, "Z_i", 'k_inf', 'k_s', 'hbar_w_ivs_h_M', 'D_ivs_h_M'...
        , "hbar_w_odp_h_M", "D_odp_h_M", "hbar_w_ivs_e_M", "hbar_w_odp_e_M" ...
        , "D_odp_e_M", "D_ivs_e_M", "rho_mass_density", "bands_type_VB", "bands_type_CB")
    end
else

    if exist('bands_transp','var') && isempty(bands_transp) == 0
        if min(CB_idx)>min(bands_transp) || max(CB_idx)<max(bands_transp)
            disp('There is an error in the attribution of the bands indexes, likely the Deformation Potentials must be extraacted for a larger energy window')
        else
             [v,p] = intersect(CB_idx,bands_transp);
             temp = D_ivs_e_M(p,p,:) ;
             D_ivs_e_M = temp;
             temp = D_odp_e_M(:,p) ;
             D_odp_e_M = temp;
             temp = hbar_w_ivs_e_M(p,p,:) ;
             hbar_w_ivs_e_M = temp;
             temp = hbar_w_odp_e_M(:,p) ;
             hbar_w_odp_e_M = temp;
        end    
    end
    if strcmp(ADP,'yes')
        save(save_filename, "Z_i", 'k_inf', 'k_s'...
       , "hbar_w_ivs_e_M", "hbar_w_odp_e_M", "D_adp_e_M", "us_sound_vel" ...
        , "D_odp_e_M", "D_ivs_e_M", "rho_mass_density", "bands_type_CB")
    else
        save(save_filename, "Z_i", 'k_inf', 'k_s'...
       , "hbar_w_ivs_e_M", "hbar_w_odp_e_M" ...
        , "D_odp_e_M", "D_ivs_e_M", "rho_mass_density", "bands_type_CB")
    end
end
