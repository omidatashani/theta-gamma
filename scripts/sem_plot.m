function sem_plot(data_matrix, x_values,c1, c2, c3)

    % Calculate mean across the first dimension
    mean_data = mean(data_matrix, 1);
    inte = 2.5;

    % Calculate SEM across the first dimension
    num_entries = size(data_matrix, 1);
    SEM_data = std(data_matrix, 0, 1) / sqrt(num_entries);

    % Define the upper and lower bounds for the SEM
    upper_bound = mean_data + SEM_data;
    lower_bound = mean_data - SEM_data;

    % Plot mean with shaded SEM
    figure;
    fill([x_values, fliplr(x_values)], [upper_bound, fliplr(lower_bound)], [c1/inte c2/inte c3/inte], 'linestyle', 'none', FaceAlpha=0.25); % Fill the area between bounds
    hold on;
    plot(x_values, mean_data, 'LineWidth', 2,'color', [c1 c2 c3]); % Plot the mean on top
    grid on;
    hold off;
end

