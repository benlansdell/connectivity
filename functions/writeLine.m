function writeLine(fh, curr_date, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual)
    %Includes only current date
    %curr_date
    fprintf(fh, '%d', curr_date);
    %cond_secs
    for idx = 1:length(cond_secs)
        fprintf(fh, '%f,', cond_secs(idx));
    end
    %electrodeBCIHpos, 
    for idx = 1:length(electrodeBCIHpos)
        fprintf(fh, '%f,', electrodeBCIHpos(idx));
    end
    %electrodeBCIVpos, 
    for idx = 1:length(electrodeBCIVpos)
        fprintf(fh, '%f,', electrodeBCIVpos(idx));
    end
    %electrodeBCIHvel, 
    for idx = 1:length(electrodeBCIHvel)
        fprintf(fh, '%f,', electrodeBCIHvel(idx));
    end
    %electrodeBCIVvel, 
    for idx = 1:length(electrodeBCIVvel)
        fprintf(fh, '%f,', electrodeBCIVvel(idx));
    end
    %electrodeBCI2vel, 
    for idx = 1:length(electrodeBCI2vel)
        fprintf(fh, '%f,', electrodeBCI2vel(idx));
    end
    %electrodeBCIDual
    for idx = 1:length(electrodeBCIDual)
        fprintf(fh, '%f,', electrodeBCIDual(idx));
    end
    fprintf(fh, '\n');
end