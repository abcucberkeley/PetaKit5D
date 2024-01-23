function cleanXmlFile(fileName)
    % Function to read a file, remove unnecessary newlines and whitespaces between specific tags,
    % and write the cleaned content back to the same file.
    
    % Step 1: Read the file content
    fileContent = fileread(fileName);
    
    % Step 2: Clean the specific tags' content to remove newlines and whitespaces
    cleanedContent = cleanSpecificTagsInText(fileContent, {'<BasePath type="relative">', '</BasePath>', '<n5 type="relative">', '</n5>'});
    
    % Step 3: Write the cleaned content back to the same file
    fid = fopen(fileName, 'w');
    fwrite(fid, cleanedContent);
    fclose(fid);
end

function cleanedText = cleanSpecificTagsInText(text, tagsToClean)
    % Helper function to remove unnecessary newlines and whitespaces between specific tags in the given text
    
    cleanedText = text;
    
    % Process each tag to clean
    for i = 1:2:numel(tagsToClean)
        openTag = tagsToClean{i};
        closeTag = tagsToClean{i+1};
        
        % Find the starting and ending indices of the tags
        startIdx = strfind(cleanedText, openTag);
        endIdx = strfind(cleanedText, closeTag);
        
        % Iterate through the occurrences of the tags
        for j = 1:numel(startIdx)
            % Extract the text between the tags
            textBetweenTags = cleanedText(startIdx(j) + length(openTag):endIdx(j) - 1);
            % Remove newlines and whitespaces from the extracted text
            cleanedTextBetweenTags = strrep(textBetweenTags, '\n', '');
            cleanedTextBetweenTags = strrep(cleanedTextBetweenTags, '\t', '');
            cleanedTextBetweenTags = strtrim(cleanedTextBetweenTags);
            % Replace the original text between tags with the cleaned version
            cleanedText = strrep(cleanedText, textBetweenTags, cleanedTextBetweenTags);
        end
    end
end
