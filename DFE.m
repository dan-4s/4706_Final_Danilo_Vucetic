function [ DFE_output ] = DFE( response )

main_cursor = response(4);
DFE_coefs = zeros(9, 1);
for i = 1: length(DFE_coefs)
    DFE_coefs(i) = -response(4+i)/main_cursor;
end

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

