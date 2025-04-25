function progress_bar(current, total)
  percent = floor((current / total) * 100);
  bar_length = 50; % length of the progress bar
  filled_length = floor(bar_length * current / total);
  bar = ['[' repmat('=', 1, filled_length) repmat(' ', 1, bar_length - filled_length) ']'];
  printf('\rProgress: %s %d%%', bar, percent);
  fflush(stdout); % force update to command line
end