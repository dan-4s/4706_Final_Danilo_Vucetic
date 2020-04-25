fprintf('=============\nQ3\n=============\n')

channel_resp_no_equ = [0.05, -0.1, 0.05, 0.58, 0.37, 0.1, 0.05, 0, -0.1, 0.04, -0.05, 0.02, 0.02];
time = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

%precursor response equalization with FFE:
ideal_precursors = [0;0;0;1];
channel_input_mat = [0.58,0.05,-0.1,0.05;0.37,0.58,0.05,-0.1;0.1,0.37,0.58,0.05;0.05,0.1,0.37,0.58];
FFE_coefs = channel_input_mat\ideal_precursors;

norm_magnitude_FFE_coefs = sum(abs(FFE_coefs));
normalized_FFE_coefs = FFE_coefs / norm_magnitude_FFE_coefs;

fprintf('The FFE precursor Matrix, ideal response, and coefficients are printed below\n');
disp(channel_input_mat);
disp(ideal_precursors);
disp(normalized_FFE_coefs);

%get new link response post-FFE
channel_output_mat = zeros(16, 4);
for i = 4 : 16
    channel_output_mat(i-3, 1) = channel_resp_no_equ(i-3);
    channel_output_mat(i-2, 2) = channel_resp_no_equ(i-3);
    channel_output_mat(i-1, 3) = channel_resp_no_equ(i-3);
    channel_output_mat(i, 4) = channel_resp_no_equ(i-3);
end

postFFE_resp = channel_output_mat*normalized_FFE_coefs;
postFFE_resp = postFFE_resp(4:16);
figure(1);
plot(time, postFFE_resp);
title('Post FFE response');
xlabel('Time in UI');
ylabel('Voltage');
grid on;


%The DFE Response now
%   Need to construct the DFE matrix. But this isn't strictly necessary
%   since DFEs can be solved by just calculating the coefficients
%   iteratively. According to Jacobs DFE notes, we can just calculate the
%   negative of the normalized ISI at each tap location. 

% We need 9 taps to cancel all the ISI
main_cursor = postFFE_resp(4);
DFE_coefs = zeros(9, 1);
for i = 1: length(DFE_coefs)
    DFE_coefs(i) = -postFFE_resp(4+i)/main_cursor;
end

%need a buffer of 3 before to get the full response
DFE_output = zeros(1, length(postFFE_resp)+9);
DFE_input = postFFE_resp;
for i = 10:length(DFE_output)
    DFE_output(i) = DFE_input(i-9);
    for j = 1:length(DFE_coefs)
        %here we multiply all preceding 9 inputs by the coef
        DFE_output(i) = DFE_output(i) + DFE_coefs(j)*DFE_output(i-j);
    end
    
end
time2 = -3:9;
DFE_output = DFE_output(10:length(DFE_output));
figure(2);
plot(time2, DFE_output);
title('Post DFE response');
xlabel('Time in UI');
ylabel('Voltage');
grid on;

fprintf('The DFE coefficients are printed below\n');
disp(DFE_coefs);

%now for the nyquist and LF responses



