function Dk_sep_scaled = filterInformationComponents(ny_data,nx_data,nz_data,norientations,norders,denom,p_vec_guess,Dk_sep,O_scaled,B_PLS,padrange,w)
%FILTERINFORMATIONCOMPONENTS Summary of this function goes here
%   Detailed explanation goes here
    Dk_sep_scaled=padarray(zeros(ny_data,nx_data,nz_data),padrange,'both');
    DK_sep_scaled = gpuArray(Dk_sep_scaled);
    
    for jj=1:norientations
        for ii=1:norders
            shift_denom=imtranslate(denom,-[p_vec_guess(ii,2,jj),p_vec_guess(ii,1,jj),p_vec_guess(ii,3,jj)]);
            numerator=padarray((Dk_sep(:,:,:,ii,jj).*conj(O_scaled(:,:,:,ii,jj).*B_PLS(ii,jj))),padrange,'both');
            Dk_sep_scaled(:,:,:,ii,jj)=numerator./(shift_denom+w^2);
        end
        
    end
end

