% Author: Alexander Rizzuto
% Date: July 16, 2023
% This code runs the Dynare code ch4_nk_taylor_minimal.mod (minimalistic
% implementation of the NK model with Taylor rule found in Gali (2008),
% Chapters 3 and 4) for various model variants (three Taylor rules) and
% compares their IRFs to a productivity (eps_a) and a m.p. shock (ni). The
% output will be stored in three different subfolders, each named after the
% Taylor rule: taylor1 is the simple Taylor rule (Ch.3, eq.25), taylor2 is
% the aggressive Taylor rule (Ch.4, eq.8), taylor3 is forward-looking
% Taylor rule (Ch.4, eq.8). Note: rates are not annualized.

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
original_path = 'ch4_nk_taylor_minimal.mod';

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

% Specify the shocks
shocknames = M_.exo_names;
varlongnames = M_.endo_names_long;
shocklongnames = M_.exo_names_long;
shocktexnames = M_.exo_names_tex;

% Set the vertical spacing between sgtitle and plots
vertical_spacing = 0.1; % Adjust this value as needed
linestyles = {'-', '--', ':', '-.'}; % Add more linestyles if needed

% Iterate over the shocks
for s = 1:numel(shocknames)
    shockname = shocknames{s};
    shocklongname = shocklongnames{s};

    % Create a figure with subplots
    figure;

    % Iterate over the filenames and load the IRF data
    for f = 1:numel(filenames)
        load(filenames{f}, 'oo_');
        irfs = oo_.irfs;

        % Get the variable names for the current shock
        varnames = fieldnames(irfs);
        varnames = varnames(endsWith(varnames, ['_', shockname]));

        % Determine the number of variables
        num_vars = numel(varnames);

        % Determine the number of rows and columns for subplot arrangement
        num_rows = ceil(sqrt(num_vars));
        num_cols = ceil(num_vars / num_rows);

        % Create a time vector for plotting
        h = size(irfs.(varnames{1}), 2);
        time = 0:h-1;

        % Adjust the vertical spacing between sgtitle and plots
        if f == 1
            ax = gca;
            ax.Position(2) = ax.Position(2) - vertical_spacing;
        end

        % Plot the IRFs in subplots
        for i = 1:num_vars
            varname = varnames{i};
            varlongname = varlongnames{i};
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
    sgtitle(['Impulse Responses to a ', shocklongname], 'FontWeight', 'bold', 'FontSize', 16);

    % Add a legend for the scenarios
    legend(suffix_names, 'Interpreter', 'none');

    % Save the figure as PNG file
    saveas(gcf, ['figure', num2str(s), '.png']);
end


clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create LaTeX File %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the filename for the LaTeX document
latex_filename = 'results.tex';

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
for s = 1:numel(shocknames)
    shockname = shocknames{s};
    shocklongname = shocklongnames{s};
    shocktexname = shocktexnames{s};
    figure_filename = ['figure', num2str(s), '.png'];

    fprintf(fid, '\\subsection{Impulse Responses to a %s $%s$}\n', shocklongname, shocktexname);
    fprintf(fid, 'The figure below shows the impulse responses to a %s for each model variant:\n\n', shocklongname);
    fprintf(fid, '\\begin{figure}[h!]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\includegraphics[width=1.2\\textwidth]{%s}\n', figure_filename);
    fprintf(fid, '\\caption{Impulse Responses to $%s$}\n', shocktexname);
    fprintf(fid, '\\end{figure}\n\n');
    fprintf(fid, '\\clearpage\n'); % Add a page break after each shock
end

% Write the LaTeX document footer
fprintf(fid, '\\end{document}\n');

% Close the LaTeX file
fclose(fid);