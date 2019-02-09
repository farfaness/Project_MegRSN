function [ labels ] = create_random_labels( nstrings )

%CREATE_RANDOM_STRING Creates random labels for creation of virtual
%channels

%   Input :  nstrings (= number of random strings) 
%   Output : labels (=cell matrix of random labels)

labels = cell(nstrings, 1);
for n=1:nstrings
    SET = char(['a':'z' '0':'9']) ;
    NSET = length(SET) ;

    N = 20 ; % pick N numbers
    i = ceil(NSET*rand(1,N)) ; % with repeat
    labels{n} = SET(i); 
end
