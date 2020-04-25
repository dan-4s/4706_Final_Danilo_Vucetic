function [ DFE_output ] = DFE( response )
%Example code for how to solve a DFE. This code was run for question 3
%where the post cursors needed to be canceled. Hence all the coeficients are
%calculated based on the known main-cursor value.
main_cursor = response(4);

%calculate the DFE coefficients. 
DFE_coefs = zeros(9, 1);
for i = 1: length(DFE_coefs)
    DFE_coefs(i) = -response(4+i)/main_cursor;
end

%calculate the DFE output given the input response. 
DFE_output = zeros(1, length(response)+9);
DFE_input = response;
for i = 10:length(DFE_output)
    DFE_output(i) = DFE_input(i-9);
    for j = 1:length(DFE_coefs)
        %here we multiply all preceding 9 inputs by the coef
        DFE_output(i) = DFE_output(i) + DFE_coefs(j)*DFE_output(i-j);
    end
    
end
DFE_output = DFE_output(10:length(DFE_output));

end

