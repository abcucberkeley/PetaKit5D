function removeTextTags(fileName)

lines = readlines(fileName);
textonly = startsWith(strip(lines),"<Text>") & endsWith(strip(lines),"</Text>");
lines(textonly) = replace(lines(textonly),["<Text>";"</Text>"],"");
writelines(lines,fileName);

end