% subfunction
function [pseudoconduction_bands, pseudovalence_bands, Ek] = shifting_bands(Ek) % %#codegen
% shifting of the bands to zero and identify the CB and VB indexes

    Emaximalextreme = zeros(1,size(Ek,4)) ;
    Eminimalextreme = zeros(1,size(Ek,4)) ;
    for id_n = 1:size(Ek,4)
        Emaximalextreme(id_n) = max(max(max(Ek(:,:,:,id_n))));
        Eminimalextreme(id_n) = min(min(min(Ek(:,:,:,id_n))));
    end

    [~,pseudoconduction_bands] = find( abs(Emaximalextreme) > abs(Eminimalextreme) );
    [~,pseudovalence_bands] = find( abs(Emaximalextreme) < abs(Eminimalextreme) );
    
    if numel(pseudoconduction_bands)>0
        new_Fermi = min(min(min(min(Ek(:,:,:,pseudoconduction_bands(1):pseudoconduction_bands(size(pseudoconduction_bands,2)))))));
        Ek = Ek - new_Fermi; % this the conduction band, or the flipped valence band, start from zero, positive values of the EF_array are into the band and negative EF values are into the gap
    end
end