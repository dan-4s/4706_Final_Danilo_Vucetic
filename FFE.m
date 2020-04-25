function [ output ] = FFE( response )
%FFE calculates the FFE response of an input channel response
%the normalized coefficients were calculated in another piece of code that is not
%presented here. See Q3_4706_Final.m
normalized_FFE_coefs = [ 0.591311905168975; -0.145323987006326;  0.172465378697213; -0.0908987291274862];

%this code just iterates through the input response calculating the 
%output values at each time step. This can also be done with a matrix. 
%See Q2_4706_final.m for examples of this. 
input = zeros(1,length(response) + 6);
input(4:length(input)-3) = response;
output = zeros(1, length(response));
for i = 1:length(response)+3
    output(i) = input(i:i+3) * normalized_FFE_coefs;
end
output = output(4:length(output));


end

