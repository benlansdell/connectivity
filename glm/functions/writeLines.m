function writeLines(fh, startdate, enddate, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual)
    %Includes startdate but not enddate
    for curr_date = startdate:(enddate-1)
        writeLine(fh, curr_date, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual);
    end
end