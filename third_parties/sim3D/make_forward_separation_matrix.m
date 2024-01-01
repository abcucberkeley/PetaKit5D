function [sep_matrix]=make_forward_separation_matrix(nphases,norders,lattice_period,phase_step)
    %Generate the forward matrix required to separate OTF orders from an image
    %acquired under multiple phases of illumination.  This matrix will need
    %to be inverted prior to multiplying by the raw image OTFs to separate
    %the components.
    
    %After this step, separation of orders will be [0, +1, -1, +2, -2,...]
    delta_phi=phase_step/lattice_period*2*pi; %Convert phase step in nm to phase step as a fraction of 2pi
    
    sep_matrix=zeros(nphases,norders);
    
    for ii=1:nphases
            counter=0;
            phase_sign=1;
        for jj=1:norders
            sep_matrix(ii,jj)=exp(1i*phase_sign*counter*(ii-1)*delta_phi);
            counter=counter+mod(jj,2);
            phase_sign=(-1+2*mod(jj,2));
        end
    end
    
    %Reorganize so that orders will be [..-2,-1,0,+1,+2...]
    sep_matrix=[sep_matrix(:,end:-2:1),sep_matrix(:,2:2:end)];