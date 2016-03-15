function [ output ] = readNodes( file )
    if nargin == 1
        fileID = fopen(file,'r');
        formatSpec = '%d %f %f %f';
        sizeA = [4 Inf];
        A = fscanf(fileID,formatSpec,sizeA);
        output=A(2:4,:);
    end
end

