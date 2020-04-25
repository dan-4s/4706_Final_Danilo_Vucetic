fprintf('=============\nQ2c\n=============\n')
%question 2c: the 1 pre, 1 main, 2 post FFE
%the matrix of equations
A = [0.58,0.05,-0.1,0.05;0.37,0.58,0.05,-0.1;0.1,0.37,0.58,0.05;0.05,0.1,0.37,0.58];
ideal_solution = [0;1;0;0];
alphas = A\ideal_solution;

%verify
for i = 1 : length(A(1,:))
    fprintf('Solution at index %d, output is: %f\n', i, A(i, :)*alphas);
end


fprintf('=============\nQ2d\n=============\n')

%2d: ZFE, need 13 coefficients, 3 pre, 1 main, 9 post
ideal_soln_2 = [0; 0; 0; 1.0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
%ideal_soln_2 = [0; 0; 0; 1.0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
% B = [ 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0;
%       0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0, 0, 0, 0;
%       0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0, 0, 0;
%       0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0, 0;
%       0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0;
%       -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0;
%       0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0;
%       -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0;
%       0.02, -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0;
%       0.02, 0.02, -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05;
%       0, 0.02, 0.02, -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1;
%       0, 0, 0.02, 0.02, -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05;
%       0, 0, 0, 0.02, 0.02, -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58];
% 
% alphas_2d = B\ideal_soln_2;
alphas_2d = A\ideal_solution;
%verify
for i = 1 : length(A(1,:))
    fprintf('Solution at index %d, output is: %f\n', i, A(i, :)*alphas_2d);
end

%getting the normalized coefficients
norm_magnitude = sum(abs(alphas_2d));
normalized_alpha = alphas_2d / norm_magnitude;

fprintf('\nNormalized alpha for 2d is below:\n')
disp(normalized_alpha);

%low frequency solutions:
all_1s = [1, 1, 1, 1];
all_0s = [-1,-1,-1,-1];

lowf_1s = all_1s * normalized_alpha;
lowf_0s = all_0s * normalized_alpha;

fprintf('\nThe low frequency response of the ZFE is %f for all 1s, and %f for all 0s\n', lowf_1s, lowf_0s);

%high frequency solutions
%low frequency solutions:
alt_1s = [1, -1, 1, -1];
alt_0s = [-1, 1, -1, 1];

highf_1s = alt_1s * normalized_alpha;
highf_0s = alt_0s * normalized_alpha;

fprintf('\nThe high frequency response of the ZFE is %f for alternating 1s, and %f for alterating 0s\n', highf_1s, highf_0s);

%calculating the peaking
hf_gain = 20*log10(abs(highf_1s));
lf_gain = 20*log10(abs(lowf_1s));
peaking = hf_gain - lf_gain;

fprintf('For 2d, peaking is: %fdB\n', peaking);
fprintf('=============\nQ2e\n=============\n')
%2e: only precursors
ideal_soln_3 = [0;0;0;1];
channel_input = [0.58,0.05,-0.1,0.05;0.37,0.58,0.05,-0.1;0.1,0.37,0.58,0.05;0.05,0.1,0.37,0.58];
alphas_2e = channel_input\ideal_soln_3;

%verify
for i = 1 : length(channel_input(1,:))
    fprintf('Solution for 2e, at index %d, output is: %f\n', i, channel_input(i, :)*alphas_2e);
end

%getting the normalized coefficients
norm_magnitude_2e = sum(abs(alphas_2e));
normalized_alpha_2e = alphas_2e / norm_magnitude_2e;

fprintf('\nNormalized alpha for 2e is below:\n')
disp(normalized_alpha_2e);

%low frequency solutions:

lowf_1s_2e = all_1s * normalized_alpha_2e;
lowf_0s_2e = all_0s * normalized_alpha_2e;

fprintf('\nThe low frequency response of the ZFE is %f for all 1s, and %f for all 0s\n', lowf_1s_2e, lowf_0s_2e);

%high frequency solutions
%low frequency solutions:

highf_1s = alt_1s * normalized_alpha_2e;
highf_0s = alt_0s * normalized_alpha_2e;

fprintf('\nThe high frequency response of the ZFE is %f for alternating 1s, and %f for alterating 0s\n', highf_1s, highf_0s);

%calculating the peaking
hf_gain_2e = 20*log10(abs(highf_1s));
lf_gain_2e = 20*log10(abs(lowf_1s_2e));
peaking_2e = hf_gain_2e - lf_gain_2e;

fprintf('For 2e, peaking is: %fdB\n', peaking_2e);
fprintf('=============\nQ2f\n=============\n')

%2f: only post cursors
% ideal_2f = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% channel_2f = [0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0, 0;
%               0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0, 0;
%               0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0, 0;
%               0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0, 0;
%               0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0, 0;
%               -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05, 0;
%               0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1, 0.05;
%               -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05, -0.1;
%               0.02,  -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58, 0.05;
%               0.02, 0.02,  -0.05, 0.04, -0.1, 0, 0.05, 0.1, 0.37, 0.58];
% alphas_2f = channel_2f\ideal_2f;

%2f: postcursor 4 tap
channel_2f = [0.58, 0.05, -0.1, 0.05;
              0.37, 0.58, 0.05, -0.1;
              0.1, 0.37, 0.58, 0.05;
              0.05, 0.1, 0.37, 0.58];
ideal_2f = [1;0;0;0];
alphas_2f = channel_2f\ideal_2f;

%verify
for i = 1 : length(channel_2f(1,:))
    fprintf('Solution for 2f, at index %d, output is: %f\n', i, channel_2f(i, :)*alphas_2f);
end

%getting the normalized coefficients
norm_magnitude_2f = sum(abs(alphas_2f));
normalized_alpha_2f = alphas_2f / norm_magnitude_2f;

fprintf('\nNormalized alpha for 2f is below:\n')
disp(normalized_alpha_2f);

%low frequency solutions:

lowf_1s_2f = all_1s * normalized_alpha_2f;
lowf_0s_2f = all_0s * normalized_alpha_2f;

fprintf('\nThe low frequency response of the ZFE is %f for all 1s, and %f for all 0s\n', lowf_1s_2f, lowf_0s_2f);

%high frequency solutions

highf_1s = alt_1s * normalized_alpha_2f;
highf_0s = alt_0s * normalized_alpha_2f;

fprintf('\nThe high frequency response of the ZFE is %f for alternating 1s, and %f for alterating 0s\n', highf_1s, highf_0s);

%calculating the peaking
hf_gain_2f = 20*log10(abs(highf_1s));
lf_gain_2f = 20*log10(abs(lowf_0s_2f));
peaking_2f = hf_gain_2f - lf_gain_2f;

fprintf('For 2f, peaking is: %fdB\n', peaking_2f);
fprintf('=============\nQ2g\n=============\n')

%2g: plotting the output waveform of 1001001, cancelling out ISI of middle
%1

%the waveform will be the superposition of the single bit response given in
%the question. So we'll need to show the time from -3UI, to 9UI + 6UI which
%is 19UI in total. I will do this by using separate vectors for each input
%at a certain time, then add the times to get the superimposed response

time = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
bit1 = [0.05, -0.1, 0.05, 0.58, 0.37, 0.1, 0.05, 0, -0.1, 0.04, -0.05, 0.02, 0.02, 0, 0, 0, 0, 0, 0];
bit2 = [0, 0, 0, 0.05, -0.1, 0.05, 0.58, 0.37, 0.1, 0.05, 0, -0.1, 0.04, -0.05, 0.02, 0.02, 0, 0, 0];
bit3 = [0, 0, 0, 0, 0, 0, 0.05, -0.1, 0.05, 0.58, 0.37, 0.1, 0.05, 0, -0.1, 0.04, -0.05, 0.02, 0.02];
overall_response = bit1+bit2+bit3;


%this one is from the excel sheet I made on drive
%overall_response = [0.05,-0.15,0.1,0.68,-0.41,-0.75,0.26,-0.51,-0.95,0.3,0.23,-0.02,0.24,-0.03,-0.09,0.09,-0.09,0,0.02];
%overall_response = [0.1	-0.1	0.1	0.15	-1.41	-1.8	-0.26	-0.71	-0.85	0.36	-0.19	-1.03	-0.8	-0.59	-0.26	0.14	-0.03	0.11	0.01];
%plot of superposition
figure(1);
plot(time, overall_response);
title('Channel response by superposition of a sequence of bits: 1001001');
xlabel('Time in UI');
ylabel('Voltage');
grid on;

%I just took 2 pre and 1 post as a kinda hail mary.. but idk what the
%actual config is supposed to be like. 
channel_resp_2g = [0.26, -0.75, -0.41, 0.68;
                   -0.51, 0.26, -0.75, -0.41;
                   -0.95, -0.51, 0.26, -0.75;
                   0.3, -0.95, -0.51, 0.26];

%the ideal response needs to have -1 in it since the channel is seeing -1s
%and 1s.
ideal_resp_2g = [-1;-1;1;-1];
alphas_2g = channel_resp_2g\ideal_resp_2g;

%verify
for i = 1 : length(channel_resp_2g(1,:))
    fprintf('Solution for 2g, at index %d, output is: %f\n', i, channel_resp_2g(i, :)*alphas_2g);
end

%getting the normalized coefficients
norm_magnitude_2g = sum(abs(alphas_2g));
normalized_alpha_2g = alphas_2g / norm_magnitude_2g;

fprintf('\nNormalized alpha for 2g is below:\n')
disp(normalized_alpha_2g);

%low frequency solutions:

lowf_1s_2g = all_1s * normalized_alpha_2g;
lowf_0s_2g = all_0s * normalized_alpha_2g;

fprintf('\nThe low frequency response of the ZFE is %f for all 1s, and %f for all 0s\n', lowf_1s_2g, lowf_0s_2g);

%high frequency solutions

highf_1s = alt_1s * normalized_alpha_2g;
highf_0s = alt_0s * normalized_alpha_2g;

fprintf('\nThe high frequency response of the ZFE is %f for alternating 1s, and %f for alterating 0s\n', highf_1s, highf_0s);

%calculating the peaking
hf_gain_2g = 20*log10(abs(highf_1s));
lf_gain_2g = 20*log10(abs(lowf_0s_2g));
peaking_2g = hf_gain_2g - lf_gain_2g;

fprintf('For 2g, peaking is: %fdB\n', peaking_2g);

fprintf('=============\nQ2h\n=============\n')

%if I slide the alpha as a sliding-window over the response, I'll get what
%it should actually look like:

time_equ = [-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10];

disp('???? idk how to do this question in matlab\n');

fprintf('=============\nQ2i\n=============\n')

%assuming 3 taps is sufficient, then using the MMSE method, we can calulate
%the least-error coefficients for the FFE. 

%have to construct the matrix now:
channel_pulse_mat = zeros(23, 5);
for i = 5 : 23
    channel_pulse_mat(i-4, 1) = overall_response(i-4);
    channel_pulse_mat(i-3, 2) = overall_response(i-4);
    channel_pulse_mat(i-2, 3) = overall_response(i-4);
    channel_pulse_mat(i-1, 4) = overall_response(i-4);
    channel_pulse_mat(i  , 5) = overall_response(i-4);
end

%3 symbols were sent in this time, at -3UI, 0UI, and 3 UI
%ideal_resp_2i = [-1;-1;-1;1;-1;-1;1;-1;-1;1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
%ideal_resp_2i = [0;0;0;1;-1;-1;1;-1;-1;1;0;0;0;0;0;0;0;0;0;0;0;0;0]; %I think this is correct.. since the '0' is actually -1
ideal_resp_2i = [0;0;0;1;0;0;1;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0]; %RZ experiment
%Then using linear least squares regression (MMSE), we can calculate the
%weights:

weights = ( transpose(channel_pulse_mat) * channel_pulse_mat ) \ transpose(channel_pulse_mat) * ideal_resp_2i ;

%verify
for i = 1 : length(channel_pulse_mat(:,1))
    %fprintf('Solution for 2i, at index %d, output is: %f\n', i, channel_pulse_mat(i, :)*weights);
    fprintf('%f\n', channel_pulse_mat(i, :)*weights);
end

%getting the normalized coefficients
norm_magnitude_2i = sum(abs(weights));
normalized_weights = weights / norm_magnitude_2i;

fprintf('\nNormalized weights for 2i is below:\n')
disp(normalized_weights);

%low frequency solutions:
all_1s_2i = [1, 1, 1, 1, 1];
all_0s_2i = [-1,-1,-1,-1,-1];

lowf_1s_2i = all_1s_2i * normalized_weights;
lowf_0s_2i = all_0s_2i * normalized_weights;

fprintf('\nThe low frequency response of the ZFE is %f for all 1s, and %f for all 0s\n', lowf_1s_2i, lowf_0s_2i);

%high frequency solutions
alt_1s_2i = [1, -1, 1,-1, 1];
alt_0s_2i = [-1, 1, -1, 1,-1];

highf_1s = alt_1s_2i * normalized_weights;
highf_0s = alt_0s_2i * normalized_weights;

fprintf('\nThe high frequency response of the ZFE is %f for alternating 1s, and %f for alterating 0s\n', highf_1s, highf_0s);

%calculating the peaking
hf_gain_2i = 20*log10(abs(highf_1s));
lf_gain_2i = 20*log10(abs(lowf_0s_2i));
peaking_2i = hf_gain_2i - lf_gain_2i;

fprintf('For 2i, peaking is: %fdB\n', peaking_2i);


