function calculate_rip()
    % Define the list of files
    files = {'dummy.csv', 'dummy1.csv', 'dummy2.csv'};
    
    % Create a map to store average RIP values per file
    ripMap = containers.Map();
    
    % Loop over each file
    for i = 1:length(files)
        total_rip = 0; % Initialize the total RIP value
        count = 0; % Initialize a counter to track the number of comparisons
        
        % Compare the current file with subsequent files
        for j = i+1:length(files)
            % Process the current file
            [y_n, gq1, tq1] = correctIQ_imbalance_t(files{i}, 0, 0);
            % Process the subsequent file
            [y_n, gq2, tq2] = correctIQ_imbalance_t(files{j}, 0, 0);
            % Calculate the RIP value and accumulate
            total_rip = total_rip + rip(gq1, tq1, gq2, tq2);
            count = count + 1; 
        end
        % Calculate the average RIP value if comparisons were made
        if count > 0
            avg_rip = total_rip / count;
        else
            avg_rip = 0; % Default to 0 if no comparisons were made
        end
        
        % Store the average RIP value in the map with the filename as the key
        if avg_rip > 0
            ripMap(files{i}) = 10*log10(avg_rip);
        else 
            ripMap(files{i}) = -70;
        end
        
        % Optionally, display the average RIP value
        fprintf('Average RIP for %s: %f\n', files{i}, ripMap(files{i}));
    end
end
