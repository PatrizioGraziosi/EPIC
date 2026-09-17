% this script, as well as the bxsf_to..., must be copied, together the .sh
% file, in the modulations directory where the mod-001... dirs are.
% in case of ADP, the ADP dirs must be in the modulations dir as well.

force_metallic = 'yes' ; % to ensure the system is treated as metallic even if a small gap is detected, useful fo rnarrow gap materials or semimetals

modes_info
%
alat = 0.1;
reduce_flag = 'n'; % 'y' if SOC and the system is not magnetic
interpolation_factor = 0;

% following, together with files names, come from the launching .sh file; 
% when I want wl make it also for QE, also alat must come from the launch 
% .sh file
%E_band_window = 0.5; % in eV, range beyong the band edge where the bands are considered, it is related to what you need in transport simulations 
%labels_q = {'G', 'A', 'K', 'H', 'M', 'L' } ;


% ======================================================================= %
% The calculation bandwidth modulations for inelastic processes start at
% line 25, for the ADP start at line 460, then the data are saved and the
% tables for inelastic processes are written, around line
% ======================================================================= %

% ----------------------------------------------------------------------- %
%        This first part computes the quantities for ODP and IVS          %
% ----------------------------------------------------------------------- %

temp = importdata('../freq.dat'); % readmatrix('../freq.dat') *4.136;  % ('../freq.txt') *4.136; % read the frequencies from THz to eV
%temp(isnan(temp)) = 0;
%freq_array = nonzeros(temp);
freq_array = temp;
%disp temp

AA = readmatrix('A_matrix.txt');
a1 = AA(1,:); a2 = AA(2,:); a3 = AA(3,:);


I_c_temp = readmatrix('coord_initial.txt');

for i = size(I_c_temp,1):-1:1
I_c(i,:) = I_c_temp(i,1)*a1 + I_c_temp(i,2) * a2 + I_c_temp(i,3)*a3 ;
end


shifting_flag = 'no';
[Ek_i_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
shifting_flag = 'yes';
[Ek_i_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;


% definition of the eventual ref. value
Ref_i = max(max(max(Ek_i_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_i_unshifted(:,:,:,1))));
% E_threshold = 0.45; % eV

% definition of VB and CB indexes
[CB_idx_temp, VB_idx_temp, ~ ] = shifting_bands(Ek_i_shifted) ; % the Ek with CB edge set to zero does not interest now

for id_band = (max(size(CB_idx_temp))):-1:1
    E_temp(id_band) = min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band)))));
end
CB_edge = min(E_temp);
clear E_temp;
CB_idx = 0*CB_idx_temp;
for id_band = (max(size(CB_idx_temp))):-1:1
    if min(min(min(Ek_i_shifted(:,:,:,CB_idx_temp(id_band))))) - CB_edge < E_band_window
        CB_idx(id_band) = CB_idx_temp(id_band);
    end
end
CB_idx = nonzeros(CB_idx);

for id_band = (max(size(VB_idx_temp))):-1:1
    E_temp(id_band) = max(max(max(Ek_i_shifted(:,:,:,VB_idx_temp(id_band)))));
end
VB_edge = max(E_temp);
clear E_temp;
VB_idx = 0*VB_idx_temp;
for id_band = (max(size(VB_idx_temp))):-1:1
    if abs(max(max(max(Ek_i_shifted(:,:,:,VB_idx_temp(id_band))))) - VB_edge) < E_band_window
        VB_idx(id_band) = VB_idx_temp(id_band);
    end
end
VB_idx = nonzeros(VB_idx);

if CB_edge > VB_edge
    semiconductor = 'yes';
    Ek_i_shifted_CB = Ek_i_shifted - CB_edge;
    Ek_i_shifted_VB = Ek_i_shifted - VB_edge;
    if strcmp(force_metallic,'yes')
        semiconductor='no';
        Ek_i_shifted_CB = Ek_i_shifted;
        CB_idx = [VB_idx ; CB_idx];
    end
else
    semiconductor = 'no';
    Ek_i_shifted_CB = Ek_i_shifted;
    CB_idx = [VB_idx ; CB_idx];
end

% find minima values and positions

% initialization
CB_i = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_i(id_band).E_ave_unshifted = 0;
    CB_i(id_band).bw = 0;
end
% evaluation
for id_band = max(size(CB_idx)):-1:1

    Ek_temp = Ek_i_shifted_CB(:,:,:,CB_idx(id_band)) ;
    Ek_temp_unshifted = Ek_i_unshifted(:,:,:,CB_idx(id_band)) ;
        
    CB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp_unshifted)));
    CB_i(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp))));

end

if strcmp(semiconductor,'yes')
    %initialization
    VB_i = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_i(id_band).bw = 0;
        VB_i(id_band).E_ave_unshifted_cut = 0;
    end
    % evaluation
    for id_band = max(size(VB_idx)):-1:1
    
        Ek_temp = Ek_i_shifted_VB(:,:,:,VB_idx(id_band)) ;
        Ek_temp_unshifted = Ek_i_unshifted(:,:,:,VB_idx(id_band)) ;
        
        VB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp_unshifted)));
        VB_i(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp))));
    
    end
end

% initialization of the DP structures
DP_bw = struct();
DP_unshifted = struct();


for im = 1:n_modes
    for iq = 1:n_q
        DP_bw(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_unshifted(im,iq).CB = zeros(max(size(CB_idx)), max(size(CB_idx)));
        DP_bw(im,iq).freq = 0;
        DP_unshifted(im,iq).freq = 0;
        if strcmp(semiconductor,'yes')
            DP_bw(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
            DP_unshifted(im,iq).VB = zeros(max(size(VB_idx)), max(size(VB_idx)));
        end
        DP_bw(im,iq).d0 = 0;
        DP_unshifted(im,iq).d0 = 0;
    end
end

counter = 0;
for iq = 1:n_q
    for im = 1:n_modes
        counter = counter + 1;
        if counter <= 9
             path_to_dir = ['mod-00',num2str(counter),'/'];
        elseif counter <= 99
            path_to_dir = ['mod-0',num2str(counter),'/'];
        else
            path_to_dir = ['mod-',num2str(counter),'/'];
        end
             
        cd(path_to_dir) % the bxsf_to... files must be copied in every dir.
        shifting_flag = 'no';
        [Ek_disp_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
        shifting_flag = 'yes';
        [Ek_disp_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ; %#ok<*ASGLU>
        
        Ref_disp = max(max(max(Ek_disp_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_disp_unshifted(:,:,:,1))));


        M_c_temp = readmatrix('coord.txt'); % read the modulated coord
        for i = size(M_c_temp,1):-1:1
            M_c(i,:) = M_c_temp(i,1)*a1 + M_c_temp(i,2) * a2 + M_c_temp(i,3)*a3 ; 
        end
        
%        D_c = M_c - I_c;
%        displ = sqrt( D_c(:,1).^2 + D_c(:,2).^2 + D_c(:,3).^2 );
	displ =  sqrt( sum( (M_c-I_c).^2,2) );
        N = displ > 10*mean(displ);
       % disp(displ(N)')
	displ(N) = max(max(AA)) - displ(N) ; % 1-displ(N);
	d0 = mean(displ) ; % average displacement in Angstrom


        % initialize the structure(s)
        CB_disp = struct();
        for id_band = (max(size(CB_idx))):-1:1
            CB_disp(id_band).E_ave_unshifted = 0;
            CB_disp(id_band).bw = 0;
        end

        for id_band = (max(size(CB_idx_temp))):-1:1
            E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx_temp(id_band)))));
        end
        CB_edge_temp = min(E_temp);
        clear E_temp;
        for id_band = (max(size(VB_idx_temp))):-1:1
            E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx_temp(id_band)))));
        end
        VB_edge_temp = max(E_temp);
        if strcmp(semiconductor,'yes')
            Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
            Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

            VB_disp = struct();
            for id_band = (max(size(VB_idx))):-1:1
                VB_disp(id_band).E_ave_unshifted = 0;
                VB_disp(id_band).bw = 0;
            end

        else
            Ek_disp_shifted_CB = Ek_disp_shifted;
        end


        % the part below, until line 214, must be repeated for VB
        for id_band = max(size(CB_idx)):-1:1
            Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
            Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,CB_idx(id_band)) ;
            
            
            CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
            CB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

        end


        % compute the deformation potentials DP
        for id_band = max(size(CB_idx)):-1:1

            DP_bw(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).bw - CB_i(id_band).bw) / d0 ; % matrixes nb x nb

            DP_unshifted(im,iq).CB(id_band,id_band) = ...
                abs( CB_disp(id_band).E_ave_unshifted - Ref_disp - ...
                ( CB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 


            for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split

                if ne(id_band,id_band_2)

                    DP_bw(im,iq).CB(id_band,id_band_2) = ...
                    abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
                    - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                    DP_bw(im,iq).CB(id_band_2,id_band) = DP_bw(im,iq).CB(id_band,id_band_2);
    
                    DP_unshifted(im,iq).CB(id_band,id_band_2) = DP_bw(im,iq).CB(id_band,id_band_2);
                    DP_unshifted(im,iq).CB(id_band_2,id_band) = DP_unshifted(im,iq).CB(id_band,id_band_2);

                end

            end

        end


        % as above for VB, from line 178
        if strcmp(semiconductor,'yes')

            for id_band = max(size(VB_idx)):-1:1

                Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;
                Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,VB_idx(id_band)) ;

                VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
                VB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

             end


            % compute the deformation potentials DP
            for id_band = max(size(VB_idx)):-1:1
    
                DP_bw(im,iq).VB(id_band,id_band) = ...
                    abs( VB_disp(id_band).bw - VB_i(id_band).bw) / d0 ; % matrixes nb x nb
    
                DP_unshifted(im,iq).CB(id_band,id_band) = ...
                    abs( VB_disp(id_band).E_ave_unshifted - Ref_disp - ...
                    ( VB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 
    
    
                for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
                    if ne(id_band,id_band_2)

                        DP_bw(im,iq).VB(id_band,id_band_2) = ...
                        abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                        - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                        DP_bw(im,iq).VB(id_band_2,id_band) = DP_bw(im,iq).VB(id_band,id_band_2);
        
                        DP_unshifted(im,iq).VB(id_band,id_band_2) = DP_bw(im,iq).VB(id_band,id_band_2);
                        DP_unshifted(im,iq).VB(id_band_2,id_band) = DP_unshifted(im,iq).VB(id_band,id_band_2);

                    end
    
                end

            end        
        end


        save mode_point_specific_data
        cd ../
        
%        whos freq_array
        f_counter = (iq-1)*total_modes + im + modes_to_skip ;
        DP_bw(im,iq).freq = freq_array(f_counter);
        DP_unshifted(im,iq).freq = freq_array(f_counter);
        DP_bw(im,iq).d0 = d0;
        DP_unshifted(im,iq).d0 = d0;

    end
end



% for the mod-all case
cd mod-all/ % the bxsf_to... files must be copied in every dir.
shifting_flag = 'no';
[Ek_disp_unshifted, ~ ] = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;
shifting_flag = 'yes';
[Ek_disp_shifted, Fermi]  = bxsf_to_ELECTRA_shifting_option(fileName,material_name,alat,reduce_flag,interpolation_factor,shifting_flag) ;

Ref_disp = max(max(max(Ek_disp_unshifted(:,:,:,1)))); % mean(mean(mean(Ek_disp_unshifted(:,:,:,1))));

M_c = readmatrix('coord.txt'); % read the modulated coord
%D_c = M_c - I_c;
%displ = sqrt( D_c(:,1).^2 + D_c(:,2).^2 + D_c(:,3).^2 );
displ =  sqrt( sum( (M_c-I_c).^2,2) );
N = displ > 10*mean(displ);
% disp(displ(N)')
displ(N) = max(max(AA)) - displ(N) ; % 1-displ(N);
d0 = mean(displ) ; % average displacement in Angstrom


% initialize the structure(s)
CB_disp = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_disp(id_band).E_ave_unshifted = 0;
    CB_disp(id_band).bw = 0;
end

for id_band = (max(size(CB_idx_temp))):-1:1
    E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx_temp(id_band)))));
end
CB_edge_temp = min(E_temp);
clear E_temp;
for id_band = (max(size(VB_idx_temp))):-1:1
    E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx_temp(id_band)))));
end
VB_edge_temp = max(E_temp);
if strcmp(semiconductor,'yes')
    Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
    Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

    VB_disp = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_disp(id_band).E_ave_unshifted = 0;
        VB_disp(id_band).bw = 0;
    end

else
    Ek_disp_shifted_CB = Ek_disp_shifted;
end


% the part below, until line 214, must be repeated for VB
for id_band = max(size(CB_idx)):-1:1
    Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
    Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,CB_idx(id_band)) ;
    
    
    CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
    CB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

end


% compute the deformation potentials DP
for id_band = max(size(CB_idx)):-1:1

    all_DP_bw.CB(id_band,id_band) = ...
        abs( CB_disp(id_band).bw - CB_i(id_band).bw) / d0 ; % matrixes nb x nb

    all_DP_unshifted.CB(id_band,id_band) = ...
        abs( CB_disp(id_band).E_ave_unshifted - Ref_disp - ...
        ( CB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ; 


    for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split

        if ne(id_band,id_band_2)

            all_DP_bw.CB(id_band,id_band_2) = ...
            abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
            - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
             all_DP_bw.CB(id_band_2,id_band) = all_DP_bw.CB(id_band,id_band_2);

            all_DP_unshifted.CB(id_band,id_band_2) = all_DP_bw.CB(id_band,id_band_2);
            all_DP_unshifted.CB(id_band_2,id_band) = all_DP_bw.CB(id_band,id_band_2);

        end

    end

end


% as above for VB, from line 178
if strcmp(semiconductor,'yes')

    for id_band = max(size(VB_idx)):-1:1

        Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;
        Ek_temp_unshifted = Ek_disp_unshifted(:,:,:,VB_idx(id_band)) ;

        VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp_unshifted)));
        VB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp))));

     end


    % compute the deformation potentials DP
    for id_band = max(size(VB_idx)):-1:1

         all_DP_bw.VB(id_band,id_band) = ...
            abs( VB_disp(id_band).bw - VB_i(id_band).bw) / d0 ; % matrixes nb x nb

         all_DP_unshifted.VB(id_band,id_band) = ...
            abs( VB_disp(id_band).E_ave_unshifted - Ref_disp - ...
            ( VB_i(id_band).E_ave_unshifted - Ref_i ) ) / d0 ;


        for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
            if ne(id_band,id_band_2)

                all_DP_bw.VB(id_band,id_band_2) = ...
                abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ d0 ;
                all_DP_bw.VB(id_band_2,id_band) = all_DP_bw.VB(id_band,id_band_2);

                all_DP_unshifted.VB(id_band,id_band_2) = all_DP_bw.VB(id_band,id_band_2);
                all_DP_unshifted.VB(id_band_2,id_band) = all_DP_unshifted.VB(id_band,id_band_2);

            end

        end

    end
    all_DP_bw.d0 = d0;
    all_DP_unshifted.d0 = d0;

end
cd ../


% ----------------------------------------------------------------------- %
%                This part computes the quantities for ADP                %
% ----------------------------------------------------------------------- %

base_name = material_name;
% ! grep 'volume' OUTCAR_base | head -n 1 > temp_volume.txt ;
! grep 'volume' ADP/base/OUTCAR | head -n 1 > temp_volume.txt ;
temp = importdata('temp_volume.txt');
base_volume = temp.data ;

% ! grep 'volume' OUTCAR_plus | head -n 1 > temp_volume.txt ;
! grep 'volume' ADP/plus/OUTCAR | head -n 1 > temp_volume.txt ;
temp = importdata('temp_volume.txt');
plus_volume = temp.data ;

% ! grep 'volume' OUTCAR_minus | head -n 1 > temp_volume.txt ;
! grep 'volume' ADP/minus/OUTCAR | head -n 1 > temp_volume.txt ;
temp = importdata('temp_volume.txt');
minus_volume = temp.data ; 





% initial
% fileName = [base_name,'_base.bxsf'] ;
% ! cp ADP/base/*bxsf ./$fileName
% ! fileName2=[base_name,'_base.bxsf'] ; cp ADP/base/*bxsf $fileName2 ; echo $fileName2
! cp ADP/base/*bxsf base_temp.bxsf 
%fileName = [base_name,'_base.bxsf'] ;

shifting_flag = 'yes';
[Ek_i_shifted_1, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,[base_name,'_base'],alat,'n',0,shifting_flag) ;
[Ek_i_shifted, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,'base_temp',alat,'n',0,shifting_flag) ;
%sum(sum(sum(Ek_i_shifted-Ek_i_shifted_1)))

% % definition of VB and CB indexes
% for id_band = (max(size(CB_idx))):-1:1
%     E_temp(id_band) = min(min(min(Ek_i_shifted(:,:,:,CB_idx(id_band)))));
% end
% CB_edge = min(E_temp); 
% clear E_temp;
% 
% for id_band = (max(size(VB_idx))):-1:1
%     E_temp(id_band) = max(max(max(Ek_i_shifted(:,:,:,VB_idx(id_band)))));
% end
% VB_edge = max(E_temp); 
% clear E_temp;
% 
% if CB_edge > VB_edge
%     semiconductor = 'yes';
%     Ek_i_shifted_CB = Ek_i_shifted - CB_edge;
%     Ek_i_shifted_VB = Ek_i_shifted - VB_edge;
% else
%     semiconductor = 'no';
%     Ek_i_shifted_CB = Ek_i_shifted;
% end

% find minima values and positions

% initialization
CB_i = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_i(id_band).E_ave_unshifted = 0;
    CB_i(id_band).bw = 0;
end
% evaluation
for id_band = max(size(CB_idx)):-1:1

    Ek_temp = Ek_i_shifted_CB(:,:,:,CB_idx(id_band)) ;
      
    CB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp)));
    CB_i(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp))));

end

if strcmp(semiconductor,'yes')
    %initialization
    VB_i = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_i(id_band).bw = 0;
        VB_i(id_band).E_ave_unshifted_cut = 0;
    end
    % evaluation
    for id_band = max(size(VB_idx)):-1:1
    
        Ek_temp = Ek_i_shifted_VB(:,:,:,VB_idx(id_band)) ;
        
        VB_i(id_band).E_ave_unshifted = mean(mean(mean(Ek_temp)));
        VB_i(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp))));
    
    end
end




% plus
fileName = [base_name,'_plus.bxsf'] ;
! cp ADP/plus/*bxsf plus_temp.bxsf 
dV = abs( (plus_volume-base_volume)/base_volume );

shifting_flag = 'yes';
% [Ek_disp_shifted, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,[base_name,'_plus'],alat,'n',0,shifting_flag) ;
[Ek_disp_shifted, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,'plus_temp',alat,'n',0,shifting_flag) ;


% initialize the structure(s)
CB_disp = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_disp(id_band).E_ave_unshifted = 0;
    CB_disp(id_band).bw = 0;
end

for id_band = (max(size(CB_idx))):-1:1
    E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx(id_band)))));
end
CB_edge_temp = min(E_temp);
clear E_temp;
if strcmp(semiconductor,'yes')
    for id_band = (max(size(VB_idx))):-1:1
        E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx(id_band)))));
    end
    VB_edge_temp = max(E_temp);
end
if strcmp(semiconductor,'yes')
    Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
    Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

    VB_disp = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_disp(id_band).E_ave_unshifted = 0;
        VB_disp(id_band).bw = 0;
    end

else
    Ek_disp_shifted_CB = Ek_disp_shifted;
end

for id_band = max(size(CB_idx)):-1:1
    Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
  
    CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp)));
    CB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

end


for id_band = max(size(CB_idx)):-1:1
    D_adp_CB_plus(id_band,id_band) = ...
        abs( CB_disp(id_band).E_ave_unshifted - CB_i(id_band).E_ave_unshifted) / (1/3*dV);
    for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split
        if ne(id_band,id_band_2)
            D_adp_CB_plus(id_band,id_band_2) = ...
                abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
                    - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ ...
            (1/3*dV);
        end
    end
end

if strcmp(semiconductor,'yes')
    for id_band = max(size(VB_idx)):-1:1
        Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;

        VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp)));
        VB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

    end

    for id_band = max(size(VB_idx)):-1:1
        D_adp_VB_plus(id_band,id_band) = ...
            abs( VB_disp(id_band).E_ave_unshifted - VB_i(id_band).E_ave_unshifted) /  (1/3*dV);
        for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
            if ne(id_band,id_band_2)
                D_adp_VB_plus(id_band,id_band_2) = ...
                    abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                        - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ ...
                (1/3*dV);
            end
        end
    end
end




% minus
% fileName = [base_name,'_minus.bxsf'] ;
! cp ADP/minus/*bxsf minus_temp.bxsf 
dV = abs( (minus_volume-base_volume)/base_volume );

shifting_flag = 'yes';
% [Ek_disp_shifted, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,[base_name,'_minus'],alat,'n',0,shifting_flag) ;
[Ek_disp_shifted, ~]  = bxsf_to_ELECTRA_shifting_option(fileName,'minus_temp',alat,'n',0,shifting_flag) ;


% initialize the structure(s)
CB_disp = struct();
for id_band = (max(size(CB_idx))):-1:1
    CB_disp(id_band).E_ave_unshifted = 0;
    CB_disp(id_band).bw = 0;
end

for id_band = (max(size(CB_idx))):-1:1
    E_temp(id_band) = min(min(min(Ek_disp_shifted(:,:,:,CB_idx(id_band)))));
end
CB_edge_temp = min(E_temp);
clear E_temp;
for id_band = (max(size(VB_idx))):-1:1
    E_temp(id_band) = max(max(max(Ek_disp_shifted(:,:,:,VB_idx(id_band)))));
end
VB_edge_temp = max(E_temp);
if strcmp(semiconductor,'yes')
    Ek_disp_shifted_CB = Ek_disp_shifted - CB_edge_temp;
    Ek_disp_shifted_VB = Ek_disp_shifted - VB_edge_temp;

    VB_disp = struct();
    for id_band = (max(size(VB_idx))):-1:1
        VB_disp(id_band).E_ave_unshifted = 0;
        VB_disp(id_band).bw = 0;
    end

else
    Ek_disp_shifted_CB = Ek_disp_shifted;
end

for id_band = max(size(CB_idx)):-1:1
    Ek_temp = Ek_disp_shifted_CB(:,:,:,CB_idx(id_band)) ;
   
    
    CB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp)));
    CB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

end



for id_band = max(size(CB_idx)):-1:1
    D_adp_CB_minus(id_band,id_band) = ...
        abs( CB_disp(id_band).E_ave_unshifted - CB_i(id_band).E_ave_unshifted) / (1/3*dV);
    for id_band_2 = max(size(CB_idx)):-1:1 % inter-band evaluated from Davydov split
        if ne(id_band,id_band_2)
            D_adp_CB_minus(id_band,id_band_2) = ...
                abs( abs( CB_disp(id_band).E_ave_unshifted - CB_disp(id_band_2).E_ave_unshifted ) ...
                    - abs( CB_i(id_band).E_ave_unshifted - CB_i(id_band_2).E_ave_unshifted ) )/ ...
            (1/3*dV);
        end
    end
end
D_adp_CB = 0.5*D_adp_CB_plus+0.5*D_adp_CB_minus;




if strcmp(semiconductor,'yes')
    for id_band = max(size(VB_idx)):-1:1
        Ek_temp = Ek_disp_shifted_VB(:,:,:,VB_idx(id_band)) ;

        VB_disp(id_band).E_ave_unshifted =  mean(mean(mean(Ek_temp)));
        VB_disp(id_band).bw = abs(max(max(max(Ek_temp)))-min(min(min(Ek_temp)))); 

    end

    for id_band = max(size(VB_idx)):-1:1
        D_adp_VB_minus(id_band,id_band) = ...
            abs( VB_disp(id_band).E_ave_unshifted - VB_i(id_band).E_ave_unshifted) /  (1/3*dV);
        for id_band_2 = max(size(VB_idx)):-1:1 % inter-band evaluated from Davydov split
            if ne(id_band,id_band_2)
                D_adp_VB_minus(id_band,id_band_2) = ...
                    abs( abs( VB_disp(id_band).E_ave_unshifted - VB_disp(id_band_2).E_ave_unshifted ) ...
                        - abs( VB_i(id_band).E_ave_unshifted - VB_i(id_band_2).E_ave_unshifted ) )/ ...
                (1/3*dV);
            end
        end
    end
    D_adp_VB = 0.5*D_adp_VB_plus+0.5*D_adp_VB_minus;

end



disp(D_adp_CB_plus)
disp(D_adp_CB_minus)
disp('')
disp(D_adp_CB)
disp('')


if strcmp(semiconductor,'yes')
    disp('')
    disp(D_adp_VB_plus)
    disp(D_adp_VB_minus)
    disp('')
    disp(D_adp_VB)

    save([base_name,'_DeformationPotentials_shifted_data.mat'], "DP_bw", "DP_unshifted", "all_DP_bw", "all_DP_unshifted", "semiconductor", "CB_idx", "VB_idx", "D_adp_CB", "D_adp_VB", "CB_disp", "CB_i", "VB_disp", "VB_i")
else
    save([base_name,'_DeformationPotentials_shifted_data.mat'], "DP_bw", "DP_unshifted", "all_DP_bw", "all_DP_unshifted", "semiconductor", "CB_idx", "D_adp_CB", "CB_disp", "CB_i")
end


n_CB = max(size(CB_idx));
n_CB_inter = (n_CB^2-n_CB)/2;
n_VB = max(size(VB_idx));
n_VB_inter = (n_VB^2-n_VB)/2;
for i_b = 1:n_CB % intra-band CB
    for im = n_modes:-1:1
        for iq = n_q:-1:1         
            matrix_f(im,iq) = DP_bw(im,iq).freq ;
            matrix(im,iq) = DP_bw(im,iq).CB(i_b,i_b) ;
        end
    end
    table_csv = array2table( matrix_f, 'VariableNames',  labels_q ) ;
    csv_name = 'freqs.csv';
    writetable( table_csv, csv_name ) ;
    table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
    csv_name = ['bw_mod_CB',i_b,'.csv'] ;
    writetable( table_csv, csv_name ) ;
end

if strcmp(semiconductor,'yes')
    for i_b = 1:n_VB % intra-band VB
        for im = n_modes:-1:1
            for iq = n_q:-1:1         
                matrix(im,iq) = DP_bw(im,iq).VB(i_b,i_b) ;
            end
        end
        table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
        csv_name = ['bw_mod_VB',i_b,'.csv'] ;
        writetable( table_csv, csv_name ) ;
    end
end

for i_b = 1:n_CB % inter-band CB
    for i_b2 = 1:n_CB
        if ne (i_b,i_b)
            for im = n_modes:-1:1
                for iq = n_q:-1:1         
                    matrix(im,iq) = DP_bw(im,iq).CB(i_b,i_b2) ;
                end
            end
            table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
            csv_name = ['DavSplit_mod_CB',i_b,'_',i_b2,'.csv'] ;
            writetable( table_csv, csv_name ) ;
        end
    end
end

if strcmp(semiconductor,'yes')
    for i_b = 1:n_VB % inter-band VB
        for i_b2 = 1:n_VB
            if ne (i_b,i_b)
                for im = n_modes:-1:1
                    for iq = n_q:-1:1         
                        matrix(im,iq) = DP_bw(im,iq).VB(i_b,i_b2) ;
                    end
                end
                table_csv = array2table( matrix, 'VariableNames',  labels_q ) ;
                csv_name = ['DavSplit_mod_VB',i_b,'_',i_b2,'.csv'] ;
                writetable( table_csv, csv_name ) ;
            end
        end
    end
end



% % subfunction
% function [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands(Ek) % %#codegen
% % shifting of the bands to zero and identify the CB and VB indexes
% 
%     Emaximalextreme = zeros(1,size(Ek,4)) ;
%     Eminimalextreme = zeros(1,size(Ek,4)) ;
%     for id_n = 1:size(Ek,4)
%         Emaximalextreme(id_n) = max(max(max(Ek(:,:,:,id_n))));
%         Eminimalextreme(id_n) = min(min(min(Ek(:,:,:,id_n))));
%     end
% 
%     [~,pseudoconduction_bands] = find( abs(Emaximalextreme) > abs(Eminimalextreme) );
%     [~,pseudovalence_bands] = find( abs(Emaximalextreme) < abs(Eminimalextreme) );
%     
%     if numel(pseudoconduction_bands)>0
%         new_Fermi = min(min(min(min(Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)))))));
%         Ek = Ek - new_Fermi; % this the conduction band, or the flipped valence band, start from zero, positive values of the EF_array are into the band and negative EF values are into the gap
%     end
% end
