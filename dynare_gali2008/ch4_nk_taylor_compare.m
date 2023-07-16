% Author: Alexander Rizzuto
% Date: July 16, 2023
% This code runs the Dynare codes ch4_nk_taylor_minimal.mod (minimalistic
% implementation) or ch4_nk_taylor.mod for various variants (three Taylor 
% rules) of the NK model with Taylor rule in Chapters 3/4 of Gali (2008),
% compares their IRFs to a technology (eps_a) and a m.p. shock (ni), and 
% performs global stability analysis. A summary of the results will be 
% stored in results.tex, the remainder in three different subfolders named 
% after the Taylor rules: taylor1 is the simple Taylor rule (Ch.3, eq.25), 
% taylor2 the "aggressive" one (Ch.4, eq.8), taylor3 the forward-looking
% one (Ch.4, eq.8).
% Note: In line 26 one can decide whether to run ch4_nk_taylor.mod or
% ch4_nk_taylor_minimal.mod.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Housekeeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the name of the mod file, including the extension
original_path = 'ch4_nk_taylor.mod'; 
%original_path = 'ch4_nk_taylor_minimal.mod';

% Model variants - Set the suffix of the output files appended to the
% original mod file name 1xN cell (N being the number of model variants)
suffix_names = {'taylor1','taylor2','taylor3'};

% Model variants - Specify the equations of the model variants
equations = {'i=phi_pi*pi+phi_y*ygap+ni;',
    'i=r_n+phi_pi*pi+phi_y*ygap+ni;',
    'i=r_n+phi_pi*pi(+1)+phi_y*ygap+ni;'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Implementing User Input %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if the length of equations and suffix_names is the same
if numel(equations) ~= numel(suffix_names)
    error('Error: The number of equations and suffix names must be the same.');
end

% Remove the file extension from the path
[~, original_name, ~] = fileparts(original_path);

% Save number of model variants
n_compare = length(suffix_names);

% Copy original mod file
old_path = 'burnin.mod';
copyfile(original_path,old_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Run Dynare for Model Variants %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:length(suffix_names)

    % Set new file name with suffix
    new_path = [original_name, '_', suffix_names{idx}, '.mod'];
    movefile(old_path, new_path);

    % Write equation to file
    fid = fopen('taylor_rules.txt', 'w');
    fprintf(fid, '%s \n', equations{idx});
    fclose(fid);

    % Run Dynare for the current model variant
    dynare(new_path, 'noclearall', 'nolog');

    % Update the file name
    old_path = new_path;

end

% Revert to original file name
movefile(old_path, original_path);

% Close all graphs
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Compare IRFs of Model Variants %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specify the filenames of the .mat files
filenames = cell(1, n_compare);
for idx = 1:n_compare
    filenames{idx} = [original_name, '_', suffix_names{idx}, '/Output/', original_name, '_', suffix_names{idx}, '_results.mat'];
end

% Combine the shock names and keep only those for which there are irfs
shock_names_n = numel(M_.exo_names);
shock_names = cell(shock_names_n, 3);
shock_names(:, 1) = M_.exo_names;
shock_names(:, 2) = M_.exo_names_long;
shock_names(:, 3) = M_.exo_names_tex;
shock_names_irfs  = unique(extractAfter(fieldnames(oo_.irfs), '_eps_'));
idx_common = ismember(extractAfter(shock_names(:, 1),'eps_'), shock_names_irfs);
shock_names = shock_names(idx_common, :);


% Combine the endo names
endo_names_n = numel(M_.endo_names);
endo_names = cell(endo_names_n, 3);
endo_names(:, 1) = M_.endo_names;
endo_names(:, 2) = M_.endo_names_long;
endo_names(:, 3) = M_.endo_names_tex;

% Set the vertical spacing between sgtitle and plots
vertical_spacing = 0.1; % Adjust this value as needed
linestyles = {'-', '--', ':', '-.'}; % Add more linestyles if needed

% Iterate over the shocks
for s = 1:size(shock_names,1)
    shock_name = shock_names{s,1};
    shock_longname = shock_names{s,2};

    % Create a figure with subplots
    figure;

    % Iterate over the filenames and load the IRF data
    for f = 1:numel(filenames)
        load(filenames{f}, 'oo_');
        irfs = oo_.irfs;

        % Get the variable names for the current shock
        endo_names_irfs = fieldnames(irfs);
        endo_names_irfs = endo_names_irfs(endsWith(endo_names_irfs, ['_', shock_name]));

        % Determine the number of variables for which there are irfs
        num_endo_irfs = size(endo_names_irfs,1);

        % Get the long names for the endo variables for which there are
        % irfs
        for i = 1:num_endo_irfs
            idx = strcmp(endo_names(:, 1), extractBefore(endo_names_irfs{i,1}, '_eps_'));
            endo_names_irfs{i, 2} = endo_names{idx, 2};
        end



        % Determine the number of rows and columns for subplot arrangement
        num_rows = ceil(sqrt(num_endo_irfs));
        num_cols = ceil(num_endo_irfs / num_rows);

        % Create a time vector for plotting
        h = size(irfs.(endo_names_irfs{1}), 2);
        time = 0:h-1;

        % Adjust the vertical spacing between sgtitle and plots
        if f == 1
            ax = gca;
            ax.Position(2) = ax.Position(2) - vertical_spacing;
        end

        % Plot the IRFs in subplots
        for i = 1:num_endo_irfs
            
            varname = endo_names_irfs{i,1};
            varlongname = endo_names_irfs{i,2};

            irf = irfs.(varname);

            subplot(num_rows, num_cols, i);
            hold on;

            % Plot the IRF from the current .mat file
            plot(time, irf(:, 1:h), 'LineWidth', 2, 'LineStyle', linestyles{f}); % Use irf(:, 1:h) to ensure consistent x-axis range
            xlim([0, h-1]); % Set the x-axis limits explicitly

            % Plot titles
            title(varlongname, 'Interpreter', 'none');

            hold off;
        end
    end

    % Add a title to the figure based on the current shock
    sgtitle(['Impulse Responses to a ', shock_longname], 'FontWeight', 'bold', 'FontSize', 16);

    % Add a legend for the scenarios
    legend(suffix_names, 'Interpreter', 'none');

    % Save the figure as PNG file
    print(['figure_', original_name, shock_name], '-dpng', '-r300');
end


clc;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create LaTeX File %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the filename for the LaTeX document
latex_filename = ['results_',original_name,'.tex'];

% Open the LaTeX file for writing
fid = fopen(latex_filename, 'w');

% Write the LaTeX document header
fprintf(fid, '\\documentclass{article}\n');
fprintf(fid, '\\usepackage{graphicx}\n');
fprintf(fid, '\\usepackage{breqn}\n');
fprintf(fid, '\\begin{document}\n');

% Start the Model Variants section
fprintf(fid, '\\section{Model Variants}\n');

% Write the Taylor rules and associated figures
for idx = 1:numel(suffix_names)
    fprintf(fid, '\\subsection{Taylor Rule %s}\n', suffix_names{idx});

    % Describe the model variant
    fprintf(fid, 'This model variant (%s) consists of the following set of equations :\n\n', suffix_names{idx});

    % Load the equations from dynamic_content.tex file
    equations_filename = fullfile([original_name, '_', suffix_names{idx}], 'latex', 'dynamic_content.tex');
    equations_content = fileread(equations_filename);

    % Write the equations to the LaTeX file
    fprintf(fid, '%s', equations_content);

    % Add the GSA figure to the LaTeX file 
    figure_filename = fullfile([original_name, '_', suffix_names{idx}, '/gsa'], [original_name, '_', suffix_names{idx}, '_prior_indeterm.pdf']);
    figure_filename = strrep(figure_filename, '\', '/'); % Replace backslashes with forward slashes
    fprintf(fid, '\\begin{figure}[h!]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\includegraphics[width=1.5\\textwidth]{%s}\n', figure_filename);
    fprintf(fid, '\\caption{Global Stability Analysis (Determinacy) for %s}\n', suffix_names{idx});
    fprintf(fid, '\\end{figure}\n');
    fprintf(fid, '\\clearpage\n'); % Add a page break after each Taylor rule



end

% Start the Impulse Responses section
fprintf(fid, '\\section{Impulse Responses}\n');

% Write the IRF figures for each shock
for s = 1:size(shock_names,1)
    shock_name = shock_names{s,1};
    shock_longname = shock_names{s,2};
    shock_texname = shock_names{s,3};
    figure_filename = ['figure_', original_name, shock_name, '.png'];

    fprintf(fid, '\\subsection{Impulse Responses to a %s $%s$}\n', shock_longname, shock_texname);
    fprintf(fid, 'The figure below shows the impulse responses to a %s for each model variant:\n\n', shock_longname);
    fprintf(fid, '\\begin{figure}[h!]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\includegraphics[width=1.2\\textwidth]{%s}\n', figure_filename);
    fprintf(fid, '\\caption{Impulse Responses to $%s$}\n', shock_texname);
    fprintf(fid, '\\end{figure}\n\n');
    fprintf(fid, '\\clearpage\n'); % Add a page break after each shock
end

% Write the LaTeX document footer
fprintf(fid, '\\end{document}\n');

% Close the LaTeX file
fclose(fid);