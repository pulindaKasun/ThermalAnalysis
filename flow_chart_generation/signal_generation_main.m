%  THE MATLAB PROGRAMME FOR THE SIGNAL GENERATION TASK IN PAGE 22 IN THE
%  THESIS
%  WARNING: MENTION THE DIRECTORIES THAT YOU WANT TO SAVE THE GRAPHS ON YOUR COMPUTER

chr = ['What model function do you want to use?' newline ...
    '1.) 1st Order' newline ...
    '2.) 2nd Order' newline ...
    '3.) A2' newline ...
    '4.) P2' newline ...
    '5.) R2' newline ...
    '6.) D1' newline ...
    '7.) Exit' newline ...
    'Enter your choice: '];

model_function = "";
invalid_input = true;

while invalid_input
    input_func = input(chr);
    
    if input_func == 1
        model_function = @first_order;
        invalid_input = false;
    elseif input_func == 2
        model_function = @second_order;
        invalid_input = false;
    elseif input_func == 3
        model_function = @A2;
        invalid_input = false;
    elseif input_func == 4
        model_function = @P2;
        invalid_input = false;
    elseif input_func == 5
        model_function = @R2;
        invalid_input = false;
    elseif input_func == 6
        model_function = @D1;
        invalid_input = false;
    elseif input_func == 7
        return;
    else
        disp("Invalid Input, Please enter a value between 1 and 7 inclusive");
    end
end
    
%  Constants
R = 8.3144621;
beta = 1/6; %  The heating rate in K/sec

%  Pre-exponential factor specifications
A_initial = 10e5;
A_final = 10e14;
A = A_initial;
A_step = 10^0.2;

%  Activation energy specifications
E_initial = 100000;
E_final = 230000;
E_step = 5000;
E_array = E_initial:E_step:E_final;
length_E_array = length(E_array);

%  Temperature specifications
T_initial = 300;  %  As it is measurable 
T_final = 2500;
T_step = 0.5;  %  0.01K in the thesis
T_array = T_initial:T_step:T_final;
length_T_array = length(T_array);

DTG_matrix = zeros(length_E_array,length_T_array,2000);  %  row: Activation energy, col: Temperature, page: pre-exponential factor 
i = 1;  %  To number the pages in the 3D array
while A <= A_final
    for j=1:length_E_array
       alpha = 0;
       alpha_array = zeros(1,length_T_array); %  To draw the alpha VS T graph
       for k=1:length_T_array
           alpha_array(k) = alpha;
           diff = (A/beta)*exp(-E_array(j)/(R*T_array(k)))*model_function(alpha);
           DTG_matrix(j,k,i) = diff;
           d_alpha = diff * T_step;
           alpha = alpha + d_alpha;
       end
       
       %  TGA curve
       plot(T_array,alpha_array);
       %str_1 = 'A: ' + A + 'E: ' + E_array(j);
       title("TGA Graph        " + 'A: ' + A + ' E: ' + E_array(j));
       xlabel("Temperature");
       ylabel("The extent of the reaction");
       saveas(gcf, "E:\Research@MSE\Theoritical\Matlab_scripts\TGA\" + i + "_" + j + "_" + k + "_TGA" + ".jpg" );
       
       %  DTG curve
       plot(T_array,DTG_matrix(j,:,i));
       %str_2 = 'A: ' + A + 'E: ' + E_array(j);
       title("DTG Graph        " + 'A: ' + A + ' E: ' + E_array(j));
       xlabel("Temperature");
       ylabel("d(alpha)/dT");
       saveas(gcf, "E:\Research@MSE\Theoritical\Matlab_scripts\DTG\" + i + "_" + j + "_" + k + "_DTG" + ".jpg" );
    end
    A = A*A_step;
    i = i + 1;
end
 
%  Saving the DTG data and the related temperature and alpha data
%  IMPORTANT: MENTION THE DIRECTORIES THAT YOU WANT TO SAVE THE GRAPHS ON YOUR COMPUTER 
% save (file location) (The variable) -> Change the file location when run on your computer
save E:\Research@MSE\Theoritical\Matlab_scripts\MAT_Files\flow_chart_signal_generation\alpha_data.mat alpha_array
save E:\Research@MSE\Theoritical\Matlab_scripts\MAT_Files\flow_chart_signal_generation\Temperature_data.mat T_array
save E:\Research@MSE\Theoritical\Matlab_scripts\MAT_Files\flow_chart_signal_generation\DTG_data.mat DTG_matrix

%  Functions for the reaction models
function value = first_order(alpha)
    value = 1 - alpha;
end

function value = second_order(alpha)
    value = (1 - alpha)^2;
end

function value = A2(alpha)
    value = 2*(1 - alpha)*(reallog(1 - alpha))^0.5;
end

function value = P2(alpha)
    value = 2*alpha^0.5;
end

function value = R2(alpha)
    value = 2*(1 - alpha)^0.5;
end

function value = D1(alpha)
    value = 1/(2*alpha);
end
