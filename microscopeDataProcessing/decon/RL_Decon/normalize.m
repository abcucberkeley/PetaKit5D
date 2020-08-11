function output = normalize(inArray, scale_factor, threshold)

% infile = [];
% if ischar(inArray)
%     infile = inArray;
%     inArray=loadtiff(infile);
% end
% 
% if ischar(scale_factor)
%     scale_factor = str2double(scale_factor);
% end
% 
% if ischar(threshold)
%     threshold = str2double(threshold);
% end

temp=sort(inArray(:));
temp=inArray./temp(round(length(temp)*threshold));
temp=image_resize(temp,round(size(temp)/scale_factor));
output = uint16(temp.*2^16);

% if ~isempty(infile)
%     thumnail_folder = 'downsampled_data/';
%     
%     inds=strfind(infile, '/');
%     if any(inds)
%         l0=inds(length(inds));
%         l1=length(infile);
%         datafolder = infile(1:l0);
%         outfile = strcat(datafolder, thumnail_folder, infile(l0+1:l1));
%     else
%         datafolder = './'
%         outfile = strcat(thumnail_folder, infile);
%     end
%     mkdir([datafolder, thumnail_folder])
%     write3Dtiff(output, outfile);
% end

end