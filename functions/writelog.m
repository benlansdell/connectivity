function writelog(logfile, str)
	fid = fopen(logfile, 'a+');
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	fprintf(fid, [stamp ': ' str '\n']);
	fclose(fid);
end