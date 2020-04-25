function [ output ] = FFE( response )
%FFE calculates the FFE response of an input channel response
normalized_FFE_coefs = [ 0.591311905168975; -0.145323987006326;  0.172465378697213; -0.0908987291274862];

input = zeros(1,length(response) + 6);
input(4:length(input)-3) = response;
output = zeros(1, length(response));
for i = 1:length(response)+3
    output(i) = input(i:i+3) * normalized_FFE_coefs;
end
output = output(4:length(output));


end

